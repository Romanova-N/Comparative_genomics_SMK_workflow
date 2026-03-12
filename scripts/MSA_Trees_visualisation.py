#This script provides the basic code to perform visualisation of MSA and phylogenetic trees construction
#What inputs we need to run it:
#1) dataset which contains enhancer name, log-value of the tree, tree type (vertebrates, mammals, primates, mixed in our case) and monophyly check (True/False)
#2) Alignment FASTA files sorted by foder according to tree type (for MSAs)
#3) Phylogenetic trees in newick format sorted by foder according to tree type (for trees)
#4) Dataset of identified motifs within aligned sequences (for MSAs)

#If you doesn't have any of that, please, adjust the code in a way of your usage our obtained all data running all snakemake files and _Main_analysis file!!
#The code in a way it is build visualisation for all provided files

#PLEASE curate the pathways:
input_dataset='path/to/input/dataset' #FILL!

trees_path = 'path/to/trees' #FILL!
tree_plot_output_path = 'output/path' #FILL!
os.makedirs(tree_plot_output_path, exist_ok=True)

alignment_input = 'path/to/alignment/fasta/files' #FILL!
plot_msa_output_path = 'output/path' #FILL!
motifs_path='path/to/identified/motifs' #FILL!
os.makedirs(plot_msa_output_path, exist_ok=True)

#input dataset
monophyly_check = pd.read_csv(input_dataset, sep=';', header=0)

#libraries
import os
import pandas as pd
import glob
import re
import gzip
from Bio import Phylo
import numpy as np
from copy import deepcopy
from collections import defaultdict
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from pymsaviz import MsaViz
from Bio import AlignIO, Align

#functions
def read_fasta(path):
    with open(path) as fasta:
        fasta = fasta.read().strip().split('>')[1:]
        fasta = [f.split('\n', 1) for f in fasta]
        fasta = {f[0]: f[1].replace('\n','').strip() for f in fasta}
    return fasta

def get_domesticated_positions(aligned_fasta):
    dom_sp = ['Felis', 'Bos', 'Canis', 'Equus', 'Sus', 'Ovis']
    other_sp = ['Homo', 'Mus', 'Pan', 'Gorilla', 'Pongo', 'Oryctolagus', 'Rattus', 'Bison']
    
    domesticated_positions = []

    for i in range(len(aligned_fasta['Homo'])):
        dom_nucls = [aligned_fasta[sp][i] for sp in dom_sp]
        other_nucls = [aligned_fasta[sp][i] for sp in other_sp]
        
        if (
            len(set(dom_nucls)) <= 2 
            and len(set(other_nucls)) <= 2
            and len(set(dom_nucls).intersection(set(other_nucls))) == 0
            ):
            domesticated_positions.append(i+1)
    return domesticated_positions

def seq_makeup(pathway):
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio import SeqIO
    import os
    
    basename = os.path.split(pathway)
    out_dir = basename[0]
    base = basename[1].split('.')[0]
    
    aligned_fasta = pathway
    file = read_fasta(aligned_fasta)
    
    right_order = ['Homo', 'Felis', 'Canis', 'Bos', 'Ovis', 'Sus', 'Equus',  'Pan', 'Gorilla',  'Pongo', 'Mus',  'Rattus', 'Oryctolagus', 'Bison'] #should be changed to your sample!!!
    reordered_file = {k: file[k].upper() for k in right_order}
    
    back_to_fasta = []
    for k, v in reordered_file.items():
        sequence = Seq(str(v))
        back_to_fasta.append(SeqRecord(sequence, id=str(k), description=""))

    with open(f'{out_dir}/{base}.fasta', 'w') as handle:
        SeqIO.write(back_to_fasta, handle, 'fasta')
    
    return print('Files are successfully made up and ready to shine in your analysis')


#Phylogenetic tree plots
for row in monophyly_check.itertuples(index=False):
    enh = row.Enhancer
    log_value = row.Log_value
    folder = row.tree_type

    
    tree_file = f"{trees_path}/{enh}.contree"
    tree = Phylo.read(tree_file, "newick")

    tree.ladderize()
    clado = deepcopy(tree)
    for clade in clado.find_clades():
        if clade.branch_length is not None and clade.branch_length > 0:
            clade.branch_length = clade.branch_length ** 0.1

    phylo_labels_rev = {
        'red':  ['Pan', 'Homo', 'Pongo', 'Gorilla'],
        'blue': ['Canis', 'Sus', 'Equus', 'Bos', 'Ovis', 'Felis']
    }
    phylo_labels_col = {label: color for color, labels in phylo_labels_rev.items() for label in labels}

    fig, ax = plt.subplots(figsize=(8, 10))

    legend_handles = [
        Patch(facecolor='red',  edgecolor='black', label='primates'),
        Patch(facecolor='blue', edgecolor='black', label='domesticated'),
        Patch(facecolor='black', edgecolor='black', label='wild')
    ]
    ax.legend(handles=legend_handles, title='Coloring', loc='upper right')

    Phylo.draw(clado, do_show=False, label_colors=phylo_labels_col, axes=ax)

    ax.set_title(f"{enh} enhancer — phylogenetic tree")
    
    fig.text(
    0.5, 0.96,
    f"log-likelihood: {log_value}",
    ha='center', va='top'
    )

    plt.tight_layout()
    plt.savefig(f'{tree_plot_output_path}/{enh}.png', dpi=300)
    plt.close(fig)
    
    
#MSAs plots
for row in monophyly_check.itertuples(index=False):
    enh = row.Enhancer
    folder = row.tree_type
    motifs = pd.read_csv(f"{motifs_path}/{enh}_motifs_mapped.csv", sep=";", header=0)
    
    aligned_fasta = f"{alignment_input}/{enh}_untrimmed.fasta"
    msa = seq_makeup(aligned_fasta)
    

    mv = MsaViz(msa, wrap_length=100, show_grid=True, color_scheme="Nucleotide", show_consensus=True)
    # Adding transcription binding sites
    for row in motifs.query('not motif_name.str.startswith("ZNF")').itertuples(index=False):
        TFBS_start = row.enh_start
        TFBS_stop  = row.enh_stop
        TFBS_name = row.motif_name

        if TFBS_name.upper().startswith("FOX"):
            color = "red"
        elif TFBS_name.upper().startswith("HOX"):
            color = "blue"  
        elif TFBS_name.upper().startswith("PAX"):
            color = "green"
        else:
            color = "black"

        mv.add_text_annotation((TFBS_start, TFBS_stop), f"{TFBS_name}", text_color=color, range_color=color)
    
    #Adding domestication associated SNPs
    mv.add_markers(get_domesticated_positions(read_fasta(msa)))
    
    mv.savefig(f'{plot_msa_output_path}/{enh}.png')