# scATAC-seq Data Analysis

## Introduction

This repository covers the analysis of scATAC-seq data produced by 10X genomics sequencing. This entry uses the data available from the [ArchR](https://www.nature.com/articles/s41588-021-00790-6) publication. Specifically, we use the scATAC-seq data for the cell line mixing experiments. The SRA information is shown below:

<pre>
SRX9633387: HighLoading
SRX9633388: LowLoading
</pre>

___

## Important Resources

[10X scATAC Paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7299161/#SD2)

[How do you convert a BAM file you found on GEO to the fastq for realignment?](https://support.10xgenomics.com/docs/bamtofastq#header)

[How do you name the Input FASTQ files for 10x pipelines that you downloaded from the internet?](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/using/fastq-input)

[How do you prepare SRA data for CellRanger?](https://kb.10xgenomics.com/hc/en-us/articles/115003802691-How-do-I-prepare-Sequence-Read-Archive-SRA-data-from-NCBI-for-Cell-Ranger-)

[SRA to Seurat walkthrough](https://bioinformaticsworkbook.org/dataAnalysis/RNA-Seq/Single_Cell_RNAseq/Chromium_Cell_Ranger.html#gsc.tab=0)

[Where can I download 10X reference genome data?](https://support.10xgenomics.com/single-cell-atac/software/downloads/latest#releasenotes)

[Can I use an ARC reference with my scATAC-seq data instead of the scATAC-seq specific reference? Yes!](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/advanced/references#arc_atac)

___

## Contents

Each section broadly cover different aspects of scATAC-seq data analysis, some more reliably than others...

### I. [Downloading 10X Data](docs/data_download.md)

### II. [Alignment and Counting with CellRanger](docs/cellranger_alignment_counting.md)

### III. Quality Control

### IIII. [Pseudo-bulk Fragment Analysis](docs/fragment_analysis.md)

