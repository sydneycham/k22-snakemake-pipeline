__author__ = "Sydney Hamilton"
__license__ = "MIT"
__email__ = "hamiltsy@ohsu.edu"

"""
CellRanger read alignment
"""

import os
import pandas as pd
import csv
import re

configfile: '/home/exacloud/gscratch/NikolovaLab/hamiltsy/k22_pipeline/cite_seq/config/config.yaml'

# SCC240606ON_Gene_Expression_library_24_00330_S7_L001_R1_001.fastq.gz

# Read the TSV file
df = pd.read_csv('alignment_samplesheet.tsv', sep='\t')

# Filter the rows where TYPE contains PROTEIN and then filter out each id
protein_df = df[df['TYPE'].str.upper().str.contains('PROTEIN')]

# use multiple samples and see what it outputs here (type)
samples = config['samples']


# def filenames(wildcards):
#     # Convert wildcards.samples to integer
#     #print(type(wildcards.samples))
#     #print(protein_df.loc[protein_df['id']] == wildcards.samples)
#     sample_id = int(wildcards.samples)

#     # Filter with pandas
#     filtered_df = protein_df[protein_df['id'] == sample_id]
#     print(protein_df['id'])
    
#     if filtered_df.empty:
#         return []
    
#     r1 = filtered_df['R1'].values[0]
#     r1_filename = r1.split("/")[-1].split(".")
#     r1_filename = r1_filename[0]
#     print(r1_filename)
    
#     r2 = filtered_df['R2'].values[0]
#     r2_filename = r2.split("/")[-1].split(".")
#     r2_filename = r2_filename[0]
#     print(r2_filename[0])
    
#     return (r1_filename,r2_filename)

def reads(wildcards):
    # Convert wildcards.samples to integer
    #print(type(wildcards.samples))
    #print(protein_df.loc[protein_df['id']] == wildcards.samples)
    sample_id = int(wildcards.samples)

    # Filter by sample id with pandas
    filtered_df = protein_df[protein_df['id'] == sample_id]
    print(protein_df['id'])
    
    if filtered_df.empty:
        return []
    
    # read1 path
    r1 = filtered_df['R1'].values[0]
    print(r1)
    
    # read2 path
    r2 = filtered_df['R2'].values[0]
    print(r2)
    
    # returning our paths to use in the input for our multiplex rule
    return (r1,r2)   

rule all:
    input:
        # Dynamically generate input paths based on sample IDs
        expand("{sample}_Multiplexing_Captures/SCC240606ON_Surface_protein_library_24_00279_S8_L001_R1_001_dup.fastq.gz", sample=samples),
        expand("{sample}_Multiplexing_Captures/SCC240606ON_Surface_protein_library_24_00279_S8_L001_R2_001_dup.fastq.gz", sample=samples),
        expand("cmo{sample}.csv", sample=samples),
        expand("final_config{sample}.csv", sample=samples),
        #expand("cellranger_log/patient_{sample}.txt", sample=samples)
        expand("test_{sample}.txt", sample=samples)

# Duplicates the protein fastqs and removes the first line LOOK AT THEM AGAIN. 
rule multiplex:
    input:
        r1 = lambda wildcards: reads(wildcards)[0],
        r2 = lambda wildcards: reads(wildcards)[1] 
    output:
        read1 = "{samples}_Multiplexing_Captures/SCC240606ON_Surface_protein_library_24_00279_S8_L001_R1_001_dup.fastq.gz",
        read2 = "{samples}_Multiplexing_Captures/SCC240606ON_Surface_protein_library_24_00279_S8_L001_R2_001_dup.fastq.gz"

    shell:
        """
        mkdir -p {wildcards.samples}_Multiplexing_Captures
        zcat {input.r1} | sed '1,4d' | gzip > {output.read1}
        zcat {input.r2} | sed '1,4d' | gzip > {output.read2}
        """


# generates cell multiplexing oligo csv with sample_id
rule cmo:
    output:
        cmo_csv = "cmo{samples}.csv"
    run:
        with open(f"cmo{wildcards.samples}.csv", 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['id','name','read','pattern','sequence','feature_type'])
            for treatment in config['treatments']:
                treatment_sequence = config['treatment_sequence'][treatment]
                writer.writerow([f"{config['samples']}_{treatment}", f"{config['samples']}_{treatment}", config['read'], config['pattern'], treatment_sequence, 'Multiplexing Capture'])


rule final_config:
    input:
        read1 = "{samples}_Multiplexing_Captures/SCC240606ON_Surface_protein_library_24_00279_S8_L001_R1_001_dup.fastq.gz",
        read2 = "{samples}_Multiplexing_Captures/SCC240606ON_Surface_protein_library_24_00279_S8_L001_R2_001_dup.fastq.gz"
    output:
        final_config = "final_config{samples}.csv"
    run:
        with open(f"final_config{wildcards.samples}.csv", 'w', newline='') as f:
            writer = csv.writer(f)
            
            # Write gene-expression section
            writer.writerow(['[gene-expression]'])
            writer.writerow(['reference', config['reference_path']])
            writer.writerow(['cmo-set', f"cmo{wildcards.samples}.csv"])
            writer.writerow(['create-bam', 'true'])
            writer.writerow([])  # Empty line

            # Write feature section
            writer.writerow(['[feature]'])
            writer.writerow(['reference', config['feature_csv_path']])
            writer.writerow([])  # Empty line

            # Write libraries section
            writer.writerow(['[libraries]'])
            writer.writerow(['fastq_id', 'fastqs', 'feature_types'])
            # Can I change the fastq_id here to something else like "Gene Expression, Multi..."
            writer.writerow(['SCC240606ON_Gene_Expression_library_24_00330', config['fastq_path'], 'Gene Expression'])
            writer.writerow(['SCC240606ON_Surface_protein_library_24_00279', config['fastq_path'], 'Antibody Capture'])
            writer.writerow(['SCC240606ON_Surface_protein_library_24_00279', os.path.abspath(os.path.dirname(input.read1)), 'Multiplexing Capture'])
            writer.writerow([])  # Empty line

            # Write samples section
            writer.writerow(['[samples]'])
            writer.writerow(['sample_id', 'cmo_ids'])
            for treatment, treatment_sequence in config['treatment_sequence'].items():
                writer.writerow([treatment, f"{config['samples']}_{treatment}"])

rule cellranger:
    input:
        final_config = "final_config{sample}.csv"
    output:
        #log = "cellranger_log/patient_{sample}.txt"
        "test_{sample}.txt"
    params:
        cellranger = "/home/exacloud/gscratch/NikolovaLab/software/cellranger-8.0.0/cellranger",
        transcriptome = "/home/exacloud/gscratch/NikolovaLab/genomes/refdata-gex-GRCh38-2020-A",
        id = "{sample}",
        cores = 24,
        mem = 64
    shell:
        """
        {params.cellranger} multi \
            --id cellranger_output_{wildcards.sample} \
            --csv {input.final_config} \
            --localcores {params.cores} \
            --localmem {params.mem} > dry_run.txt 
        touch test_{wildcards.sample}.txt
        """