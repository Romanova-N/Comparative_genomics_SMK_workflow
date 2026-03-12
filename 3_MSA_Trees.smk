#The third snakemake which creates MSA and consensus trees from prepared sequences
#NOTE: This snakemake takes longer than all others, check timing in README.md

#======= TO RUN:
# snakemake -s 3_MSA_Trees.smk --use-conda --cores 32 -n -p
# snakemake -s 3_MSA_Trees.smk --use-conda --cores 32 &> logs/3_MSA_trimming_Phylo.log

## libraries
import os
import pandas as pd
from snakemake.io import glob_wildcards

configfile: 'config.yaml'

## global wildcards
tree_dir = "output/3_Phylogen_trees"
os.makedirs(tree_dir, exist_ok=True)

#if you want define specific groups

# enh_groups = pd.read_csv("data/2_MSA_input_fasta_files/_enhancer_groups.csv", sep=";", header=0)
# targets = {"full_sample", config['target_name']}
# selected_enh = enh_groups[enh_groups["tree_type"].isin(targets)]
# enhancers = selected_enh.Enhancer.astype(str).tolist()

## functions
def getMSAs_names(wildcards):
    ck_dir = checkpoints.msa_trim.get(**wildcards).output[0]
    names = [os.path.splitext(fn)[0] for fn in os.listdir(ck_dir)
        if fn.endswith(".fasta") and "_untrimmed" not in fn and os.path.isfile(os.path.join(ck_dir, fn))]
    return names


def getTrees_names(wildcards):
    names = getMSAs_names(wildcards)
    return expand(os.path.join(tree_dir, "{enhancer}_tree.contree"), enhancer=names)

## expected output
rule all:
    input:
        "output/_trees_enhancers.txt",


#================= PART 1: MSA (batch checkpoint) =================#
checkpoint msa_trim:
    """
    Run t_coffee on each enhancer FASTA and trim with trimal.
    Note: if trimal filters out 1 or more sequences, this sample will not proceed to 
    phylogenetic tree construction. The trimmed file will be deleted, but untrimmed is kept.
    """
    threads: 16
    input:
        expand("data/2_MSA_input_fasta_files/{enhancer}.fasta", enhancer=enhancers)
    output:
        directory("output/3_MSA_results")
    params:
        method = "promo_pair@EP@GOP@-60@GEP@-1,muscle_msa"
    shell:
        r"""
        outdir={output}
        untrim="$outdir/untrimmed"
        mkdir -p "$outdir" "$untrim"
        shopt -s nullglob

        # Parallelised alignment+trimming
        process_one() {{
            fasta="$1"
            [ -s "$fasta" ] || exit 0

            enh="$(basename "$fasta" .fasta)"
            out_trim="$outdir/${{enh}}.fasta"
            out_aln="$untrim/${{enh}}_untrimmed.fasta"
            out_html="$untrim/${{enh}}_untrimmed.html"

            t_coffee \
                -seq "$fasta" \
                -type=dna \
                -method {params.method} \
                -n_core 4 \
                -output fasta_aln \
                -outfile "$out_aln" \
                -outorder=input >/dev/null 2>&1 || true


            trimal -in "$out_aln" -fasta -automated1 -out "$out_trim" -htmlout "$out_html" >/dev/null 2>&1 || true

            # if at least 1 species seq was filtered out -> delete trimmed
            orig_len=$(grep -c '^>' "$out_aln" || true)
            trim_len=$(grep -c '^>' "$out_trim" || true)

            if [ "$orig_len" -ne "$trim_len" ]; then
                rm -f "$out_trim"
                exit 0
            fi

            exit 0
        }}

        export -f process_one
        export outdir untrim
        
        mv -f ./*.dnd "$untrim/" 2>/dev/null || true
        """

#================= PART 2: Trees (batch checkpoint) =================#

rule phylogenetic_trees:
    """
    Build trees for all trimmed MSAs. The rule starts only after all MSAs are generated.
    """
    threads: 2
    input:
        msa=lambda wc: os.path.join(checkpoints.msa_trim.get(**wc).output[0], f"{wc.enhancer}.fasta")
    output:
        contree=os.path.join(tree_dir, "{enhancer}_tree.contree")
    params:
        iqargs = "-bb 1000 -m GTR -mem 4G",
        prefix=os.path.join(tree_dir, "{enhancer}_tree")
    shell:
        r"""
        iqtree -s {input.msa} -pre {params.prefix} {params.iqargs} -nt {threads} || true

        [ -f {output.contree} ] || : > {output.contree}
        """

#================= PART 3: List of trees =================#

rule list_tree_enhancers:
    """
    Collecting names of enhancers for which we got phylogenetic tree
    """
    input:
        getTrees_names
    output:
        "output/_trees_enhancers.txt"
    shell:
        r"""
        mkdir -p $(dirname {output})
        : > {output}
        for t in {input}; do
            [ -s "$t" ] || continue
            base="$(basename "$t" .contree)"
            echo "$base" >> {output}
        done
        """
