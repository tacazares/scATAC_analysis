# Pseudo-bulk Fragment Analysis

This section covers analysis of pseudo-bulk level scATAC-seq data.

## Important Resources
[ArchR paper](https://www.nature.com/articles/s41588-021-00790-6)

## Fragment Files

The [scATAC-seq fragment files](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/output/fragments) are like bed files, but with a `.tsv` extension. The fragment files are representative of paired end sequencing reads.

### Example Fragment File

```tsv
chr1	10315	10531	TCAGGGCCACTAAACC-1	1
chr1	49477	49855	CACTGAATCTATCTAC-1	1
chr1	181119	181610	GTTATTCAGTCGGGAT-1	1
chr1	181250	181470	GCAGCCAGTTCCGCGA-1	1
chr1	181323	181489	CTGGGACCACTCCCAT-1	1
chr1	181348	181389	CTAGCGGAGATGTTGA-1	2
chr1	181356	181470	GAGTGAGTCGCTATAG-1	1
chr1	181370	181470	GCTTAAGCACTGCTTC-1	1
chr1	181390	181458	TTACCCGAGAACGCCA-1	1
chr1	181405	181453	CAACGTAGTAGCGGTA-1	1
```

### Column Description

| Column Number |     Name    |                                                                       Description                                                                       |
|:-------------:|:-----------:|:-------------------------------------------------------------------------------------------------------------------------------------------------------:|
| 1             | chrom       | Reference genome chromosome of fragment                                                                                                                 |
| 2             | chromStart  | Adjusted start position of fragment on chromosome.                                                                                                      |
| 3             | chromEnd    | Adjusted end position of fragment on chromosome. The end position is exclusive, so represents the position immediately following the fragment interval. |
| 4             | barcode     | The 10x cell barcode of this fragment. This corresponds to the CB tag attached to the corresponding BAM file records for this fragment.                 |
| 5             | readSupport | The total number of read pairs associated with this fragment. This includes the read pair marked unique and all duplicate read pairs.                   |

### Tn5 specific correction

The fragments have been adjusted for the 9bp overhang created by Tn5 transposition, so the ends of each fragment are representative of the Tn5 cut sites. My analysis is primarily interested in inferring the exact Tn5 site that was bound, so we will convert the paired-end sequencing fragments to Tn5 sites.

This analysis uses 10 cell types that were used in the [ArchR](https://www.nature.com/articles/s41588-021-00790-6) paper. Each cluster of cells was identified and extracted as individual files containing the fragments associated with each cell type.

## Converting Fragments to Tn5 insertion sites

I convert the scATAC-seq fragment files to bed files, where each entry represents a Tn5 cut site that has been windowed to 40bp. I chose 40bp based on the x-ray crystallography of Tn5 that shows it has ~19bp wide binding site.

```
Fragment       ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓  
Tn5 IS         █                █
Slop 20bp    █████            █████
```

Below are 3 entries from the fragments file HEK293T cells:

```bash
$ head -n 3 fragments_Granja2021_293T.bed
chr1	10315	10531	TCAGGGCCACTAAACC-1	1
chr1	49477	49855	CACTGAATCTATCTAC-1	1
chr1	181119	181610	GTTATTCAGTCGGGAT-1	1
```

You can convert the fragments to Tn5 sites by splitting the start and stop positions into separate entries and windowing around them. The last position of the fragment is exclusive, so we -1 from the stop position to find the exact base that represents the cut site. This is because [BED intervals are 0-start half-open coordinates](https://genome-blog.soe.ucsc.edu/blog/2016/12/12/the-ucsc-genome-browser-coordinate-counting-systems/), so they do not include the last position in the range when extracting sequences or values.

```bash
# Tn5 sites based on fragment start position
chr1	10315	10316	TCAGGGCCACTAAACC-1	1
chr1	49477   49478	CACTGAATCTATCTAC-1	1
chr1	181119  181120	GTTATTCAGTCGGGAT-1	1

# Tn5 sites based on fragment stop position
chr1	10530	10531	TCAGGGCCACTAAACC-1	1
chr1	49854   49855	CACTGAATCTATCTAC-1	1
chr1	181609  181610	GTTATTCAGTCGGGAT-1	1
```

The python script [fragmentsToTn5.py](../scripts/scATAC_fragments_tsv_to_tn5_bed.py) will perform the fragment splitting and can be found in the [scripts](../scripts/) directory.

I window the cut sites using `bedtools slop`:

```bash
bedtools slop -i {cut_sites} -g hg38.chrom.sizes -b 20 > windows_cut_sites.bed
```

```bash
chr1	10295	10336	TCAGGGCCACTAAACC-1
chr1	49457	49498	CACTGAATCTATCTAC-1
chr1	181099	181140	GTTATTCAGTCGGGAT-1
```

### Table of Counts

| Cell Type | Tn5 Cut Sites | Fragments | Barcode Count |
|:---------:|---------------|:---------:|---------------|
| HEK293T   | 26,299,054      | 13,149,527  | 1,756          |
| Jurkat    | 20,076,486      | 10,038,243  | 1,670          |
| HeLa      | 18,727,198      | 9,363,599   | 2,385          |
| GM12878   | 17,966,182      | 8,983,091   | 1,291          |
| K562      | 11,436,140      | 5,718,070   | 1,494          |
| THP-1     | 7,500,868       | 3,750,434   | 1,059          |
| MCF-7     | 4,346,222       | 2,173,111   | 1,208          |

## Convert Tn5 sites to coverage tracks

The next step is to convert the Tn5 cut sites into a coverage track that is scaled to the number of insertion sites.

We first need to calculate the scaling factor that will be used with `bedtools genomecov`. The calculation is broken up, so you can easily see how each number was derived.

```bash
# Count the number of fragments or entries
frag_count=$(cat ${file} | wc -l)

# Divide 1 by the fragment count
reads_factor=$(bc -l <<< "1/${frag_count}")

# Multiply by the scaling factor
rpm_factor=$(bc -l <<< "${reads_factor} * 1000000")
```

Then use `bedtools genomecov`:

```bash
bedtools genomecov -i - -g hg38.chrom.sizes -bg -scale ${rpm_factor} > ${bedgraph_file}
```

A bedgraph only shows the # of bed intervals at each position. You will lose the UMI information and any meta data when making a bedgraph or bigwig.

```bash
chr1	10295	10336	0.0380242
chr1	10510	10551	0.0380242
chr1	49457	49498	0.0380242
```

Once you have a [bedgraph](https://genome.ucsc.edu/goldenpath/help/bedgraph.html) file that contains the genome coverage information, you can compress it to a [bigwig](https://genome.ucsc.edu/goldenPath/help/bigWig.html) to save space and increase speed when loading/transferring. Use [bedGraphToBigWig](https://anaconda.org/bioconda/ucsc-bedgraphtobigwig) for quick conversion. Below is an example of how to convert a bedgraph file to bigwig.

```bash
bedGraphToBigWig {bedgraph_file} {chrom_sizes} {bigwig_file}
```

## Calling peaks on Tn5 cut sites

We can call peaks for areas of enriched chromatin accessibility using the Tn5 cut sites generated from the fragment files. We call peaks using [MACS2](https://github.com/macs3-project/MACS) and the bed file of Tn5 cut sites.

### MACS2 parameters

`--extsize 40` Slop size to extend intervals from the left most side
`-f BED` Input is BED format
`--nomodel` Use ext size above
`--keep-dup all` Keep all duplicates
`-B` Saves the bedgraph
`--SPMR` Normalizes the signal to number of reads in the file

### Example MACS2 code

```bash
macs2 callpeaks -t {tn5_bed} --extsize 40 -f BED --name {output_name} --nomodel --keep-dup all -g hs -B --SPMR
```

## Generate an indexed fragments file for input to Seurat, Sinto, Signac....

This information was directly copied from the site for [Sinto](https://timoast.github.io/sinto/scatac.html).

```bash
# sort, compress, and index
sort -k1,1 -k2,2n fragments.bed > fragments.sort.bed
bgzip -@ 8 fragments.sort.bed
tabix -p bed fragments.sort.bed.gz
```

## UMAP of scATAC-seq Tn5 cut site peaks

We can generate a union or oracle peak set from the peaks that we call individually on each cell type. 

You can combine the peaks using `bedtools`, `cat`, and `cut`.

```bash
cat ./Granja_2021_scatac/peaks/*bed | cut -f1,2,3 | bedtools sort | bedtools merge > Granja_2021_union_peaks.bed
```

### Generate Counts Matrix

We will generate the Tn5 counts matrix by intersecting the oracle peak set with each of the sample Tn5 cut sites. This will give us a count of Tn5 sites per peak for each of the cell types.

```bash
for file in /Users/caz3so/scratch/20210831_maxATAC_scATAC/data/hg38/tn5_slop20/*gz;
do
bedtools intersect -a Granja_2021_7CellTypes_union_peaks.bed -b ${file} -c > /Users/caz3so/scratch/20210831_maxATAC_scATAC/outputs/hg38/umap/counts/`basename ${file} _IS_slop20.bed.gz`_counts.bed
done
```