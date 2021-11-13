# Alignment and counting with CellRanger

## Important Resources

[How do you prepare SRA data for CellRanger?](https://kb.10xgenomics.com/hc/en-us/articles/115003802691-How-do-I-prepare-Sequence-Read-Archive-SRA-data-from-NCBI-for-Cell-Ranger-)

[SRA to Seurat walkthrough](https://bioinformaticsworkbook.org/dataAnalysis/RNA-Seq/Single_Cell_RNAseq/Chromium_Cell_Ranger.html#gsc.tab=0)

[Where can I download 10X reference genome data?](https://support.10xgenomics.com/single-cell-atac/software/downloads/latest#releasenotes)

[Can I use an ARC reference with my scATAC-seq data instead of the scATAC-seq specific reference? Yes!](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/advanced/references#arc_atac)

## Running CellRanger

The data was aligned based on the workflow described in the [10X tutorial](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/using/count).

### Example code: running CellRanger

```bash
# Load the module
module load cellrangerATAC/2.0.0

# Change to the output directory
cd /scratch/caz3so/Granja2021_scATAC/output

# Run CellRanger ATAC count

# High Loading Sample
cellranger-atac count --id=SRX9633387 \
--reference=./genome_inf/hg38/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
--fastqs=/scratch/caz3so/Granja2021_scATAC/input/fastq/CellLine_HighLoading/CellLine_High_MissingLibrary_1_H2MCLBGXC \
--sample=bamtofastq \
--localcores=60 \
--localmem=120

# Low Loading Sample
cellranger-atac count --id=SRX9633388 \
--reference=./genome_inf/hg38/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
--fastqs=/scratch/caz3so/Granja2021_scATAC/input/fastq/CellLine_HighLoading/CellLine_High_MissingLibrary_1_H2MCLBGXC \
--sample=bamtofastq \
--localcores=60 \
--localmem=120
```

### Run Statistics

#### High Loading

|          Feature         |   Memory/Time/Count  |
|:------------------------:|:--------------:|
| CPU time :               | 200476.00 sec. |
| Max Memory :             | 73786 MB       |
| Average Memory :         | 43048.44 MB    |
| Total Requested Memory : | 768000.00 MB   |
| Delta Memory :           | 694214.00 MB   |
| Max Swap :               | 2 MB           |
| Max Processes :          | 123            |
| Max Threads :            | 810            |
| Run time :               | 13041 sec.     |
| Turnaround time :        | 13053 sec.     |

#### Low Loading

|          Feature         |   Memory/Time/Count  |
|:------------------------:|:--------------:|
| CPU time :               | 203892.44 sec. |
| Max Memory :             | 64872 MB       |
| Average Memory :         | 32376.43 MB    |
| Total Requested Memory : | 256000.00 MB   |
| Delta Memory :           | 191128.00 MB   |
| Max Swap :               | 2 MB           |
| Max Processes :          | 124            |
| Max Threads :            | 856            |
| Run time :               | 14806 sec.     |
| Turnaround time :        | 14811 sec.     |