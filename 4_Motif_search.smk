#The forth and the last snakemake file which works with raw enhancer fasta files to identify motif regions either only for reference species or for all
#NOTE: By default, the workflow identifies motif regions only in a reference sequence.  If the snakemake map motifs to all species, it's quite long!
#NOTE: coordinates of output motifs are adjusted to MSA files (i.o. so biased by indels)


#======= TO RUN:
#snakemake -s 4_Motif_search.smk --use-conda -j 16 -n -p #dry-run
#snakemake -s 4_Motif_search.smk --use-conda -j 16 &> logs/4_Motif_search.log

## libraries
import os
import re
import glob
import pandas as pd
from pymemesuite.fimo import FIMO 
from pymemesuite.common import MotifFile 
from pathlib import Path
from pymsaviz import MsaViz
from Bio import AlignIO, Align
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

configfile: 'config.yaml'

## global wildcards
# list of enhancers with available trees (and MSAs), formed by rule list_tree_enhancers in the 3d snakemake file
aligned_enh = pd.read_csv('output/_trees_enhancers.txt', header=None)[0].astype(str).tolist()

JASPAR_URL = "https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt"


## functions
def read_fasta(path):
    """Read a multi-FASTA file into {first_token_of_header: sequence} dict."""
    with open(path) as fasta:
        fasta = fasta.read().strip().split('>')[1:]
        fasta = [f.split('\n', 1) for f in fasta]
        fasta = {f[0].split()[0]: f[1].replace('\n','').strip() for f in fasta}
    return fasta

def generate_mapping(raw_fasta, aligned_fasta, organism=None, name_trimmed=True):
    """Build raw->aligned coordinate mapping (1-based positions, gaps in aligned)."""
    if organism is None:
        raise ValueError("organism must be provided")
    raw_seq = str(raw_fasta[organism.split()[0] if name_trimmed else organism])

    aln_key = organism.split()[0] if name_trimmed else organism
    aln_seq = str(aligned_fasta[aln_key])

    mapping = {}
    raw_pos = 0
    for aln_pos, ch in enumerate(aln_seq, start=1):
        if ch == '-':
            continue
        raw_pos += 1
        mapping[raw_pos] = aln_pos
        if raw_pos >= len(raw_seq):
            break
    return mapping

def map_coordinates(enh_start, enh_end, raw_fasta, aligned_fasta, organism, **kwargs):
    """Map raw start/stop to aligned coordinates; return (mapped_start, mapped_stop)."""
    m = generate_mapping(raw_fasta, aligned_fasta, organism, **kwargs)
    s = m.get(int(enh_start))
    e = m.get(int(enh_end))
    return s, e


## expected output
rule all:
    input:
        expand("output/4_Motifs_remapped_for_msa_plot/{enhancer}_motifs_mapped.csv", enhancer=aligned_enh)


#================= PART 1: Prepare JASPAR database =================#

rule download_meme_library:
    """Download JASPAR (MEME format) as a reproducible step."""
    output:
        "data/4_JASPAR_motifs_database.txt"
    shell:
        r"""
        mkdir -p $(dirname {output})
        curl -L "{JASPAR_URL}" -o {output}
        """

#================= PART 2: FIMO motif search =================#

rule motif_search:
    """ 
    With available database of vertebrates motifs, the rule look through fasta sequences before multiple sequence alignment
    and identify all motifs in reference sequence (!!!). As an output, it provides .csv file per enhancer which lists all motifs in the sequences
    in which seqience it was found, statistical scoring, etc.
    """  
    input:
        motifs="data/4_JASPAR_motifs_database.txt",
        fasta="data/2_MSA_input_fasta_files/{enhancer}.fasta"
    output:
        "output/4_Motifs_per_enhancer/uncurated/{enhancer}_motifs.txt"
    params:
        ref_latin = config["ref_latin"],
        base_dir = "output/4_Motifs_per_enhancer/uncurated",
        thresh = "1e-4",
        markov_order = 1
    shell:
        r"""
        mkdir -p {params.base_dir} logs/fimo

        outdir="{params.base_dir}/{wildcards.enhancer}"
        mkdir -p "$outdir"

        ref_fa="$outdir/{wildcards.enhancer}.ref.fasta"
        bg="$outdir/{wildcards.enhancer}.bg"

        # Extract only the reference organism sequence
        awk -v ref="{params.ref_latin}" '
            BEGIN {{p=0}}
            /^>/ {{
                p = (index($0, ref) > 0)
            }}
            p {{print}}
        ' {input.fasta} > "$ref_fa"

        #  Build background model (Markov)
        # fasta-get-markov writes background file; keep it in per-enhancer outdir
        fasta-get-markov -m {params.markov_order} "$ref_fa" "$bg" >/dev/null 2>&1

        # adjusting p-value threshold (if seq length > 1kb)
        ref_len=$(
        awk '
            /^>/ {next}
            { gsub(/[ \t\r\n]/, "", $0); len += length($0) }
            END { print len+0 }
        ' "$ref_fa"
        ) #calculation of the sequence length

        p0="{params.thresh}"   # e.g. 1e-4
        p_eff=$(
        awk -v L="$ref_len" -v p="$p0" '
            BEGIN {
            if (L > 1000) printf "%.10g", (p * 1000.0 / L);
            else printf "%.10g", p;
            }
        '
        ) #recalculation of the p-value threshold

        # Run FIMO; do not fail the rule if it fails
        fimo --oc "$outdir" --thresh $p_eff --bgfile "$bg" {input.motifs} "$ref_fa" >/dev/null 2>&1

        # Pick the produced FIMO table (format may vary by version)
        src=""
        if [ -f "$outdir/fimo.tsv" ]; then
          src="$outdir/fimo.tsv"
        elif [ -f "$outdir/fimo.txt" ]; then
          src="$outdir/fimo.txt"
        fi

        # If no output produced -> create empty file, but do not crash
        if [ -z "$src" ] || [ ! -s "$src" ]; then
          echo "FIMO produced no result table for {wildcards.enhancer}"
          : > {output}
          exit 0
        fi

        mv "$src" {output}
        """

rule curating_meme_output:
    """ 
    Initial output of fimo doesn't provide normal naming of found motifs. The dataset contains only motifs IDs which
    is not really convenient in analysis proceeding. Thus, the rule takes JASPAR database, where each motif is named as
    "MOTIF MA0069.1 PAX6", dissect IDs and proper name, goes through all output csvs from the previous rule and replace
    IDs by proper name. 
    """  
    input:
        motifs_txt = "output/4_Motifs_per_enhancer/uncurated/{enhancer}_motifs.txt",
        meme_lib   = "data/4_JASPAR_motifs_database.txt"
    output:
        motifs_csv = "output/4_Motifs_per_enhancer/{enhancer}_motifs.csv"
    run:
        # dictionary "MAxxxx.x -> motif_name"
        mapping = {}
        with open(input.meme_lib, encoding="utf-8") as fh:
            for line in fh:
                m = re.match(r'^MOTIF\s+(\S+)(?:\s+(.+))?$', line.strip())
                if m:
                    motif_id = m.group(1)
                    alt = (m.group(2) or "").strip()
                    if alt:
                        mapping[motif_id] = alt

        # reading FIMO output
        df = pd.read_csv(input.motifs_txt, sep=None, engine="python")

        rename_map = {
            "#pattern name": "motif_id",
            "motif_id": "motif_id",           
            "sequence name": "organism",
            "matched sequence": "sequence",
            "p-value": "p_value",
            "q-value": "q_value",
            "start": "start",
            "stop": "stop",
        }
        df = df.rename(columns={k: v for k, v in rename_map.items() if k in df.columns})

        # column with readble motif name
        if "motif_id" in df.columns:
            df["motif_name"] = df["motif_id"].map(mapping).fillna(df["motif_id"])
        else:
            df["motif_name"] = ""

        df.to_csv(output.motifs_csv, sep=';', index=False)


#================= PART 3: Motifs coordinates for alignments =================#

rule remapping_motifs:
    """ 
    Motif search can be performed only with classic fasta files meanwhile we need to have start-end positions for aligned
    files in order to use it for visualisation (MSA plots).
    """
    input:
        raw_fasta     = "data/2_MSA_input_fasta_files/{enhancer}.fasta",
        aligned_fasta = "output/3_MSA_results/untrimmed/{enhancer}_untrimmed.fasta",
        motifs_csv    = "output/4_Motifs_per_enhancer/{enhancer}_motifs.csv"
    output:
        motifs_remapped = "output/4_Motifs_remapped_for_msa_plot/{enhancer}_motifs_mapped.csv"
    run:
        raw_fasta = read_fasta(input.raw_fasta)
        aligned_fasta = read_fasta(input.aligned_fasta)
        motifs = pd.read_csv(input.motifs_csv, sep=';', header=0)

        for col in ['organism', 'start', 'stop']:
            if col not in motifs.columns:
                raise ValueError(f"Required column '{col}' is missing in {input.motifs_csv}")

        motifs['mapped_start'] = None
        motifs['mapped_stop']  = None

        for row in motifs.itertuples():
            enh_start = getattr(row, 'start')
            enh_stop  = getattr(row, 'stop')
            organism  = getattr(row, 'organism')

            try:
                new_start, new_stop = map_coordinates(enh_start, enh_stop, raw_fasta, aligned_fasta, organism)
            except Exception:
                new_start, new_stop = (None, None)

            motifs.at[row.Index, 'mapped_start'] = new_start
            motifs.at[row.Index, 'mapped_stop']  = new_stop

        motifs.to_csv(output.motifs_remapped, sep=";", index=False)
