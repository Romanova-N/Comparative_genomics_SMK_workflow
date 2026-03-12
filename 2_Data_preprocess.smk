#This is the second snakemake file which works with data uploaded by 1_Download_data.smk file

#======= TO RUN:
#snakemake -s 2_Data_preprocess.smk -n -p  #dry-run
#snakemake -s 2_Data_preprocess.smk --use-conda --cores 1 &> logs/2_Data_preprocess.log

## libraries
import os
import glob
import gzip
from tqdm import tqdm
import pandas as pd
import math
from collections import defaultdict
import bioframe as bf
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.AlignIO import MultipleSeqAlignment


configfile: 'config.yaml'

## global wildcards
species_database = pd.read_csv(config['sample_list'], sep=';', header=0, index_col=False)
ref = config['ref_organism']
want_all_annotations = bool(config.get('annotations_needed', False))


assembly_list = species_database['Assembly'].tolist()
assembly_list_without_ref = [x for x in assembly_list if x != ref]
latin_list = species_database['Latin'].astype(str).str.replace(' ', '_', regex=False).tolist()
species_to_latin = dict(zip(assembly_list, latin_list)) #for set_sequence function (line 119)

#preparing all directories
os.makedirs("output/2_enhancers_genes_pairs", exist_ok=True)
os.makedirs("data/2_comparative_datasets", exist_ok=True)
os.makedirs("data/2_MSA_input_fasta_files", exist_ok=True)

## functions
def set_sequence(coord, sp, fasta_dict, species_to_latin=species_to_latin):
    """
    the function takes coordinate in a format"chromosome\tstart\tend", species name and a dictionary with a genome procceded 
    with SeqIO.to_dict and gives the same coordinates, sequence of these coordinates and species name as an output
    """
    if pd.isna(coord):
        return None
    chrom, start_i, end_i = coord.split('\t')
    start, end = int(start_i ) - 1, int(end_i)
    
    rec = fasta_dict.get(chrom)
    if rec is None:
        return None
    sequence = rec.seq[start:end]

    latin_name = species_to_latin[sp]
    return [coord, sequence, latin_name]

## expected output
rule all:
    input:
        expand("output/2_enhancers_genes_pairs/{assembly}.enh_gene_pairs.csv",
                    assembly=assembly_list if want_all_annotations else [ref]),
            "data/2_comparative_datasets/enhs_without_duplicates.csv",
            "data/2_comparative_datasets/enh_sequences_dataset.pkl",
            "data/2_MSA_input_fasta_files/_enhancer_groups.csv",


#============= PART 1: Enhancer-gene pairs=============#

rule prepare_coord_to_bioframe:
    """
    To measure distances between the enhancer and TSS region, we need to dissect all
    TSSs positions accounting a DNA strand. This will be used in the enhancer_gene_pairing rule
    """  
    input:
        gff_csv="data/1_GFF_gene_annotations/{assembly}.genes_coordinates.csv"
    output:
        bioframe_files="data/2_bioframe_tmp_files/{assembly}.bioframe_tmp.csv"
    run:
        df = pd.read_csv(input.gff_csv, sep=";", header=0)
        df['start'] = df['start'].where(df['strand'] == '+',other=df['end'])
        new_gff = df.get(['Chr', 'start']).copy()
        new_gff['start'] = new_gff['start'].astype(int)
        new_gff['end'] = new_gff['start'] + 1 #just to technically get a "region" as further software expects

        new_gff.to_csv(output.bioframe_files, sep=";", index=False)


rule enhancer_gene_pairing: 
    """
    The rule construct enhancer-gene pairs based ONLY on the distance between enhancer median position and gene's TSS
    """  
    input:
        bioframe_files="data/2_bioframe_tmp_files/{assembly}.bioframe_tmp.csv",
        gff_csv="data/1_GFF_gene_annotations/{assembly}.genes_coordinates.csv",
        enh_coord=lambda wc: (f"data/1_liftOver/{wc.assembly}.enh_coordinates.csv") if wc.assembly != ref else config['ref_enh']
    output:
        enh_gene_pairs="output/2_enhancers_genes_pairs/{assembly}.enh_gene_pairs.csv"
    params:
        ref = ref
    run:
        #getting the data
        genes_tss = pd.read_csv(input.bioframe_files, sep=';', header=0)
        genes_ann = pd.read_csv(input.gff_csv, sep=';', header=0)
        if wildcards.assembly == params.ref:
            enh_df = pd.read_csv(input.enh_coord, sep="\t", header=None, usecols=[0,1,2],
                         names=["chrom","start","end"]) #annotation only for the reference, DEFAULT
        else:
            enh_df = pd.read_csv(input.enh_coord, sep=";", header=0) #annotation for each species, requires all GFF in data!

        #defining the closest TSS for each enhancer
        all_closest = []
        for chrom, gsub in genes_tss.groupby('Chr'):
            esub = enh_df[enh_df['chrom'] == chrom]
            if not esub.empty:
                df1 = esub[['chrom','start','end']]
                df2 = gsub[['Chr','start','end']].rename(columns={'Chr':'chrom'})
                closest = bf.closest(df1, df2, suffixes=("_E","_G"))
                all_closest.append(closest)
        if all_closest:
            df = pd.concat(all_closest, ignore_index=True)
        else:
            df = pd.DataFrame()
        
        #Substitution of TSS coordinates by actual gene coordinates (DNA strand is accounted)
        merged_df_pos = genes_ann.merge(df[["chrom_E", "start_E", "end_E", "chrom_G", "start_G", "distance"]], 
                            left_on=["Chr", "start"], right_on=["chrom_G", "start_G"], how="inner") #for genes on + strand
        merged_df_neg = genes_ann.merge(df[["chrom_E", "start_E", "end_E", "chrom_G", "start_G", "distance"]], 
                            left_on=["Chr", "end"], right_on=["chrom_G", "start_G"], how="inner") #for genes on - strand
        
        #final dataset with gene info and related enhancer
        merged_df = pd.concat([merged_df_pos, merged_df_neg], ignore_index=True)
        result = merged_df[["seqid", "gene", "Chr", "start", "end", "strand", "start_E", "end_E", "distance"]].copy()
        result.to_csv(output.enh_gene_pairs, sep=';', index=False)


#============= PART 2: Enhancer prepared to MSA =============#
rule enhancer_comparative_dataset:
    """
        Based on reference-organism coordinates, we build a wide table: 1st column — reference enhancer; following columns — target-species enhancers.
        One row represents which homologous regions a reference enhancer has across other species. If a cell is NA,
        there was no corresponding region found in the target species.
    """   
    input:
        ref_path = config['ref_enh'],
        liftovers=expand("data/1_liftOver/{sp_name}.enh_coordinates.csv", sp_name=assembly_list_without_ref)
    output:
        enh_without_duplicates="data/2_comparative_datasets/enhs_without_duplicates.csv"
    params:
        ref = ref
    run:
        #getting coordinated of reference enhancers
        final_df = pd.read_csv(input.ref_path, sep='\t', header=None, usecols=[0,1,2], names=['ref_chrom','ref_start','ref_end'])
        final_df["ref_chrom"] = final_df["ref_chrom"].astype(str)
        final_df[["ref_start","ref_end"]] = final_df[["ref_start","ref_end"]].astype(int)

        #adding a ref species column and multiidex from ref coordinates
        final_df[params.ref] = final_df['ref_chrom'].astype(str) + "\t" + final_df['ref_start'].astype(str) + "\t" + final_df['ref_end'].astype(str)
        final_df.set_index(['ref_chrom', 'ref_start', 'ref_end'], inplace=True)

        #adding all other species coordinates ignore all cases when enhancer got mapped in several regions (drop duplicates)
        for sp in input.liftovers:
            df = pd.read_csv(sp, sep=";", header=0)
            sp_name = os.path.basename(sp).split(".enh_coordinates.csv")[0]
            df["chrom"] = df["chrom"].astype(str)
            df[["start","end"]] = df[["start","end"]].astype(int)

            df[sp_name] = (df["chrom"].astype(str) + "\t" + df["start"].astype(str) + "\t" + df["end"].astype(str))

            subdf = df.set_index(["ref_chrom", "ref_start", "ref_end"])[[sp_name]] #reindexing
            subdf = subdf[~subdf.index.duplicated(keep=False)]
            final_df = final_df.join(subdf, how="left") #merging based on multiindex

        #save
        final_df.to_csv(output.enh_without_duplicates, sep=";", index=True, header=True)

rule enhancer_sequences_dataset:
    """
    Here we take previously prepared dataset with all corresponding enhancers. Using coordinates, we find the sequence of each of them
    and put into the dataframe cell as a SeqRecord object.
    """   
    input:
        ref_gene_annotation=f"output/2_enhancers_genes_pairs/{ref}.enh_gene_pairs.csv",
        enh_without_duplicates="data/2_comparative_datasets/enhs_without_duplicates.csv",
        sp_fasta_files=expand("data/1_FASTA_genomes/{assembly}.fa.gz", assembly=assembly_list)
    output:
        enhancer_sequences_dataset="data/2_comparative_datasets/enh_sequences_dataset.pkl"
    run:
        ref_annotation=pd.read_csv(input.ref_gene_annotation, sep=';', header=0)
        enh_dataset=pd.read_csv(input.enh_without_duplicates, sep=';', header=0, index_col=['ref_chrom','ref_start','ref_end'])

        #enumerating of genes to avoid duplicates in names
        genes = list(ref_annotation['gene'])

        counts = defaultdict(int)
        unique = []
        for g in genes:
            counts[g] += 1
            suffix = f"_{counts[g]}"
            unique.append(g + suffix)   
        ref_annotation['genes'] = pd.Series(unique) #so we got SOX9_1 and SOX9_2 if 2 enhancers got annotated to SOX9 gene

        #joining gene names and enhancer coordinates per species dataset by multiindex
        annotated_enh = ref_annotation[['genes', 'Chr', 'start_E', 'end_E']].rename(columns={'Chr': "ref_chrom", 'start_E': 'ref_start', 'end_E': 'ref_end'})
        annotated_enh.set_index(['ref_chrom', 'ref_start', 'ref_end'], inplace = True)

        genes_species = annotated_enh.join(enh_dataset, how="inner")
        genes_species.set_index('genes', inplace=True)

        #for each coordinate we get sequences from genome assemblies, saved as SeqRecord per cell
        full_data = pd.DataFrame()
        for sp_fasta in tqdm(input.sp_fasta_files):
            sp_name = os.path.basename(sp_fasta).split(".fa.gz")[0]
            if sp_name not in genes_species.columns:
                continue
            with gzip.open(sp_fasta, "rt") as handle:
                fasta_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))

            full_data[sp_name] = genes_species[sp_name].apply(lambda coord: set_sequence(coord, sp_name, fasta_dict))
            del fasta_dict
        
        data = full_data.T
        data.to_pickle(output.enhancer_sequences_dataset)
        
checkpoint preparing_fasta_files_for_MSA:
    """
    Final preparation for MSA: using the generated dataset containing enhancer sequences, it generate FASTA files per enhancer
    Structure:
    > species name
    sequence
    """
    input:
        enhancer_sequences = "data/2_comparative_datasets/enh_sequences_dataset.pkl"
    params:
        fasta_dir = "data/2_MSA_input_fasta_files"   
    output:
        fasta_dir = directory("data/2_MSA_input_fasta_files"),
        groups_csv = "data/2_MSA_input_fasta_files/_enhancer_groups.csv"
    run:
        enh_sequences = pd.read_pickle(input.enhancer_sequences)

        total_number_of_species = len(enh_sequences.index)

        rows = [] #collect info about enhancer type

        #the loop makes fasta record (">" + seq) per each enhancer
        for gene in enh_sequences.columns:
            records = []
            group_check = set()

            for sp in enh_sequences.index:
                entry = enh_sequences.at[sp, gene]
                if not entry:
                    continue
                seq_obj, latin = entry[1], entry[2]
                records.append(SeqRecord(seq_obj, id=str(latin), description=""))
                group_check.add(str(latin).split(' ')[0])

            if not records:
                continue

            # write FASTA into the single folder
            fasta_path = os.path.join(output.fasta_dir, f"{gene}.fasta")
            SeqIO.write(records, fasta_path, "fasta")

            #based on species set we assign each enhancer to particular group
            #it doesn't affect fasta files or the way it stored, but the df
            #will be used for wildcards in the next rules
            target_group = set(config['target_group'])
            if len(records) == total_number_of_species:
                tree_type = "full_sample"
            elif group_check == target_group:
                tree_type = config['target_name']
            else:
                tree_type = "other"
            rows.append({"Enhancer": gene, "tree_type": tree_type})

        # save the mapping CSV
        group_dataset = pd.DataFrame(rows, columns=["Enhancer", "tree_type"])
        group_dataset.to_csv(output.groups_csv, index=False, sep=";")


