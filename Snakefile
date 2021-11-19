configfile: "./tests/data/config.yaml"

# This script was adapted from ChIP and ATAC workflows from https://github.com/crazyhottommy/

# localrules will let the rule run locally rather than submitting to cluster
localrules: all

import pandas as pd
import os
import glob

def get_samples(input_dir):
    return glob.glob(os.path.join(input_dir, "*.bed"))

samples_df = pd.read_table(config['SAMPLES_TSV']).set_index("Sample_ID", drop=False)

ALL_SAMPLES = samples_df["Sample_ID"].unique().tolist()

rule all:
    input:  ALL_BIGWIG  + ALL_FLAGSTAT + ALL_FASTQC + ALL_PEAKS

rule fragments_to_tn5:
    input: 
        os.path.join(config["input_dir"], "{sample}.tsv.gz")
    output: 
        file_path=os.path.join(config["output_dir"],"Tn5_CutSites", "{sample}.bed.gz")
    params: 
        split=config["split_ends"],
        dir_path=config["output_dir"]
    log: 
        os.path.join(config["output_dir"], "{sample}/logs/{sample}.fragments_to_tn5.txt")
    benchmark: 
        os.path.join(config["output_dir"], "{sample}/logs/{sample}.fragments_to_tn5.benchmark.txt")
    threads: 
        4
    conda: 
        "./envs/genome_tools.yaml"
    message: 
        "Convert fragments to Tn5 sites {slop_size}: {threads} threads"
    shell:
        """
        python ./scripts/fragmentsToTn5.py -i {input} -o {params.dir_path} -n {sample}
        """

rule bedtools_slop:
    input: 
        bed=os.path.join(config["output_dir"],"Tn5_CutSites", "{sample}.bed.gz")
    output: 
        os.path.join(config["output_dir"], "{sample}", "Tn5_slop", "{sample}_slop.bed.gz")
    params: 
        blacklist_bed=config["blacklist_bed"],
        slop_size=config["slop_size"],
        chrom_sizes=config["reference_chrom_sizes"]
    log: 
        os.path.join(config["output_dir"], "{sample}/logs/{sample}.bedtools_slop.txt")
    benchmark: 
        os.path.join(config["output_dir"], "{sample}/logs/{sample}.bedtools_slop.benchmark.txt")
    threads: 
        4
    conda: 
        "./envs/genome_tools.yaml"
    message: 
        "Remove blacklisted regions, slop intervals by {slop_size}, sort, and compress with pigz: {threads} threads"
    shell:
        """
        bedtools intersect -a {bed_file} -b {blacklist_bed} -v | bedtools slop -i - -g {chrom_sizes} -b {slop_size} | sort -k 1,1 -k2,2n | pigz > {output}
        """

rule bedtools_genomecov:
    input: 
        os.path.join(config["output_dir"], "{sample}", "Tn5_slop", "{sample}_slop.bed.gz")
    output: 
        temp(os.path.join(config["output_dir"], "{sample}", "Tn5_slop", "{sample}_slop.bg"))
    params: 
        chrom_sizes=config["reference_chrom_sizes"]
    log: 
        os.path.join(config["output_dir"], "{sample}/logs/{sample}.bedtools_genomecov.txt")
    benchmark: 
        os.path.join(config["output_dir"], "{sample}", "logs", "{sample}.bedtools_genomecov.benchmark.txt")
    threads: 
        4
    conda: 
        "./envs/genome_tools.yaml"
    message: 
        "Use bedtools to calculate genome coverage for {input}: {threads} threads"
    shell:
        """
        sh ./scripts/BedToRPMBedgraph.sh {input} {output_dir} {chrom_sizes}
        """

rule bedGraphToBigWig:
    input:
        bedGraph=os.path.join(config["output_dir"], "{sample}", "Tn5_slop", "{sample}_slop.bg"),
        chromsizes=config["reference_chrom_sizes"]
    output:
        os.path.join(config["output_dir"], "{sample}", "bigwig", "{sample}.bw")
    log:
        os.path.join(config["output_dir"], "{sample}", "logs", "{sample}.bed-graph_to_big-wig.txt")
    params:
        "" # optional params string
    wrapper:
        "v0.80.1/bio/ucsc/bedGraphToBigWig"