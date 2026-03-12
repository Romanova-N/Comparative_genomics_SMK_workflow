#This snakemake workflow will prepare all initial files for analysis:
    #-enhancer coordinates
    #-genomes in .fasta format
    #-gene annotations in .gff format
#NOTE: we work only with so-called active enhancers (poised enhancers are ignored)
#NOTE: rules download_chain_files and download_genomes_fasta don't work with more than 1 core, UCSC server limitations!!

#======= TO RUN:
#snakemake  -s 1_Download_data.smk -n -p    #dry-run
#snakemake -s 1_Download_data.smk --use-conda -j 1 &> logs/1_Download_data_workflow.log

## libraries
import pandas as pd
import os

configfile: 'config.yaml'

## global wildcards
species_database = pd.read_csv(config['sample_list'], sep=';', header=0, index_col=False)
ref = config['ref_organism']
# ref_list = [ref]
want_all_annotations = bool(config.get('annotations_needed', False))
assembly_list = species_database.Assembly.tolist()
assembly_list_without_ref = [x for x in assembly_list if x != ref]

if want_all_annotations: 
    assembly_to_annotation = dict(zip(species_database.Assembly, species_database.Annotation))
else:
    assembly_to_annotation = {ref: config['ref_organism_annotation']}


## expected output
rule all:
    input:
        # liftover CSV 
        expand("data/1_liftOver/{assembly}.enh_coordinates.csv", assembly=assembly_list_without_ref),
        # FASTA genome assemblies
        expand("data/1_FASTA_genomes/{assembly}.fa.gz", assembly=assembly_list),
        # genome annotation
        expand("data/1_GFF_gene_annotations/{assembly}.genes_coordinates.csv",
            assembly=assembly_list if want_all_annotations else ref)

###======== PART 1 Lift over ========###
rule download_chain_files:
    """
    Downloading all chain files of listed assemblies
    """
    output:
        chain="data/chain_files/{assembly}.over.chain.gz"
    threads: 1
    params:
        upper_name=lambda wc: wc.assembly[0].upper() + wc.assembly[1:],
    shell:
        """
        mkdir -p data/chain_files

        url="ftp://hgdownload.soe.ucsc.edu/goldenPath/{ref}/liftOver/{ref}To{params.upper_name}.over.chain.gz"
        echo "Download $url ..., {output.chain}"
        wget -O {output.chain} "$url" --timeout=1 --tries=5
        """

rule liftover:
    """
    LiftOver  enhancers to target species using chain files.
    Resulting files have less number of lines since not all reference coordinates have homologous regions in
    target species genomes (due to indels).
    """
    input:
        active_enhs=config['ref_enh'],
        sp_chain="data/chain_files/{assembly}.over.chain.gz"
    output:
        mapped="data/1_liftOver/Act_Enhancer/{assembly}.bed",
        unmapped="data/1_liftOver/Act_Enhancer/unmapped/{assembly}.txt"
    shell:
        """
        mkdir -p data/1_liftOver/Act_Enhancer/unmapped
        liftOver -minMatch=0.1 -multiple -preserveInput {input.active_enhs} \
                 {input.sp_chain} {output.mapped} {output.unmapped}
        """


rule preprocess_liftover:
    """
    Transfer all lifted over coordinates from bed files to .csv dataset
    """
    input:
        bed="data/1_liftOver/Act_Enhancer/{assembly}.bed"
    output:
        csv="data/1_liftOver/{assembly}.enh_coordinates.csv"
    run:
        if os.path.getsize(input.bed) == 0:
            open(output.csv, "w").close()
        else:
            df = pd.read_csv(input.bed, sep="\t", header=None)
            df.columns = ["chrom","start","end","Ref_Coord","strain"]
            df[["ref_chrom","ref_start","ref_end"]] = df["Ref_Coord"].str.replace(':','-').str.split('-', expand=True)
            df.drop("Ref_Coord", axis=1, inplace=True) #enh coordinates of reference are splited in 3 columns
            df["ref_start"] = df["ref_start"].astype(int) - 1
            df["ref_end"] = df["ref_end"].astype(int)
            df.to_csv(output.csv, sep=";", index=False)


###======== PART 2 Download genome annotations (GFF) ========###
rule download_gene_annotation_gff:
    """
    It download gene annotation following IDs listed in initial file (or only reference annotation)
    """
    output:
        gff="data/1_GFF_gene_annotations/{assembly}/genomic.gff",
        report="data/1_GFF_gene_annotations/{assembly}/sequence_report.jsonl"
    params:
        annotation=lambda wildcards: assembly_to_annotation[wildcards.assembly]
    shell:
        """
        mkdir -p data/1_GFF_gene_annotations/{wildcards.assembly}
        datasets download genome accession {params.annotation} --filename data/1_GFF_gene_annotations/{wildcards.assembly}/{wildcards.assembly}.zip --include gff3,seq-report

        unzip -q data/1_GFF_gene_annotations/{wildcards.assembly}/{wildcards.assembly}.zip -d data/1_GFF_gene_annotations/{wildcards.assembly}/
        mv data/1_GFF_gene_annotations/{wildcards.assembly}/ncbi_dataset/data/{params.annotation}/genomic.gff data/1_GFF_gene_annotations/{wildcards.assembly}/
        mv data/1_GFF_gene_annotations/{wildcards.assembly}/ncbi_dataset/data/{params.annotation}/sequence_report.jsonl data/1_GFF_gene_annotations/{wildcards.assembly}/

        rm -rf data/1_GFF_gene_annotations/{wildcards.assembly}/ncbi_dataset data/1_GFF_gene_annotations/{wildcards.assembly}/{wildcards.assembly}.zip
        rm -rf data/1_GFF_gene_annotations/{wildcards.assembly}/README.md data/1_GFF_gene_annotations/{wildcards.assembly}/md5sum.txt 
        """

rule preprocess_annotations:
    """
    It preprocess initial annotation data to the convenient .csv file 
    """
    input:
        gff_raw="data/1_GFF_gene_annotations/{assembly}/genomic.gff",
        json_report="data/1_GFF_gene_annotations/{assembly}/sequence_report.jsonl"
    output:
        gff_csv="data/1_GFF_gene_annotations/{assembly}.genes_coordinates.csv"
    run:
        #importing gff file with genome annotation and extracting only protein-coding annotations
        col_names = ['seqid','source','type','start','end','score','strand','phase','attributes']
        gff = pd.read_csv(input.gff_raw, sep="\t", comment="#", header=None, names=col_names)

        genes = gff.loc[gff["type"].eq("gene"), ["seqid","start","end","strand","attributes"]].copy()
        genes["gene"] = genes["attributes"].str.extract(r"(?:^|;)gene=([^;]+)")
        genes["gene_biotype"] = genes["attributes"].str.extract(r"(?:^|;)gene_biotype=([^;]+)")
        gff2 = genes.loc[genes["gene_biotype"].eq("protein_coding"), ["seqid","gene","start","end","strand"]]

        #importing of seq reports with correspondence of NCBI sequence IDs and chromosome numbers
        report = pd.read_json(input.json_report, lines=True)

        if 'ucscStyleName' in report.columns and report['ucscStyleName'].notna().sum() > 2:
            chr_col = 'ucscStyleName'
        elif 'chrName' in report.columns:
            chr_col = 'chrName'
        else:
            raise ValueError(f"{wildcards.assembly} has not appropriate file format, download manually")

        report = report[[chr_col, 'refseqAccession']].copy()
        if chr_col == 'chrName':
            report[chr_col] = "chr" + report[chr_col].astype(str)
        report.columns = ['Chr', 'seqid']

        #merging IDs with chromosome number in a final dataset, filtering out genes on unknown chromosome (not interpretable)
        gff3 = gff2.merge(report, on=['seqid'], how='left')
        final_gff = gff3.drop(gff3.index[gff3['Chr'] == 'Un'])

        #saving a new dataset with species short name
        final_gff.to_csv(output.gff_csv, sep=";", index=False)

###======== PART 3 Download genomes (FASTA) ========###

rule download_genomes_fasta: 
    """
    It downloads all target genomes fro UCSC database in accordance with chain files version
    """
    params:
        ucsc_url="ftp://hgdownload.soe.ucsc.edu/goldenPath/{assembly}/bigZips/{assembly}.fa.gz"
    output:
        genome_fasta="data/1_FASTA_genomes/{assembly}.fa.gz"
    shell:
        """
        mkdir -p data/1_FASTA_genomes
        echo "======{wildcards.assembly}======"

        wget -O {output.genome_fasta} {params.ucsc_url}
        """