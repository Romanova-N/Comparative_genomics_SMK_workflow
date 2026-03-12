#### Comparative research of neural crest regulome in light f domestication

### Overview
The pipeline was developed to facilitate the search of corresponding DNA regions across species via UCSC LiftOver tools. Pdrticulary it was used for the research aiming to understand the impact of regualatory regions involved in neural crest dvelopment in observed domestication phenomenon.

The repository contains several snakemake files which should be run separately one by one (see running instrutions in the first lines of each file). Before starting make sure you filled up the config file with all required input information.
After filling the config, you need to run 0_Input_validation.py script, to check whethereverything is correctly assigned and your working directory contains all nessecary files to run the pipeline seamlessly.

**Step-by-step instruction:**
1) Set up the conda environment (instructions below)
1) Fill up the config file, follow the instructions inside
2) Run 0_Input_validation.py script. In case everything is fine you get "[OK] Preflight checks passed." message
3) Run the snakemake files one by one, following the nubers in file names. To do it, open a file and follow the instructions in the header.

Please, always do the dry run first and then proceed to the actual process.

====================
=====conda env======
in the direction of the project run
conda env create -f envs/main.yaml

before running the pipeline, activate the environment:
conda activate name_of_env
====================

### Short descriptions of the snakemake files
- **1_Download_data.smk**
As the name says, the main purpose to get all necessary files from different databases. It downloads chain files from UCSC database and corresponfing genome assemblies in FASTA format. 
It also downloads a genome annotation either only for the reference species (defined by config file) or for all species (than all assemblies ID in RefSeq format must be addet in the sample list!!).
Finally, the rule carries out the liftover of genome coordinated and prepare the results in BED format
Important notes:
1) This rule can be run ONLY with 1 CORE. Unfotunately, if you do too much requests to USCS website at the same time they block it.

- **2_Data_preprocess.smk**
To make the data weсв  download easier to use, we need to preprocess it. First, we annotate all reference enhancers (from the ref species) by proximation methods: the nearest TSS of a gene is assigned to the enhancer.
Then the script form a dataframe containing all corresponding region coordinates across species, so index is a reference enhancer info (gene assighned+coordinates) and each column is a species coordinates.
Then, using genome assemblies from eac species, we add actual sequences of that regions to each dataframe cell
Important notes:
1) Since one gene can have several enhancers, enhancer annotation has a number coming afte the gene name. So if we have 3 enhancers for ETS1 gene, they will be annotated like ETS1_1, ETS1_2, ETS1_3. The number doesn't
reflect the actual order of these enhancers in the genome!!
2) Sometimes, due to Lift over, one reference enhancer is splited into several ones in a target genome. It generates an enormous amount of duplicates in resulting dataframe where all species are put together.
Due to that reason, the pipeline completely delete enhancer information if it cause such an effect. Thus, in resulting dataset, we have only those enhancers which were mapped unseparately to all targeted genomes. 

- **3_MSA_Tree.smk**
This is the key part of the pipeline: it performs multiple sequence alignment per enhancer and then builds a consensus tree. The MSA uses 2 methods: the penalties for regulatory regions recommended by T-Coffee developers 
and muscle iterative algorythm. Using two methods simuntaneously increase the accuracy, accorsing to the literature. Then, the alignment proceed to the trimming (it's desirable for phylogenetic trees based on Maximum Likelihood Approach). Those alighments which lost at least one species during the trimming due to low level of similarity, don't procceed to phylogenetic tree recnstruction. Finally, the pipeline build a phylogenetic tree based on trimmed alighment. It implements GTR substitution model and maximum likelihood approach, all trees are supported by bootstrap (which altogether makes the process quite slow). In the end, the pipeline generates a _trees_enhancers.txt list of all successfully generated trees.
Note: 

- **4_Motif_search.smk**
The last pipeline perform the motif search across all raw fasta files (!!) and then recalulates resulting coordinates to alighnment coordinates. Importantly: FIMO looks for motifs only in reference species sequence. According to FIMO reccommendations, we form a appropriate background with fasta-get-markov and adjust p-value threshold if the sequence is > 1000 bp (for more details read FIMO Tutorial)


Times and resouces for 4k enhancers, 22 species:
1) 1_Download_data.smk - 30 min, STRICTLY 1 core
2) 2_Data_preprocess.smk - 10 min, multiple cpus will not improve the time
3) 3_MSA_Tree.smk - 8 h, more cpus can be provided to enhance phylogenetic trees construction
4) 4_Motif_search.smk - 6 h


In /scripts folder you will find additional MSA_Trees_visualisation which need to be adjusted for your purposes