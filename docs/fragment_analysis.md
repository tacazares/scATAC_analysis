# Pseudo-bulk Fragment Analysis

This section covers analysis of pseudo-bulk level scATAC-seq data.

## Important Resources
[ArchR paper](https://www.nature.com/articles/s41588-021-00790-6)

## Fragment Files

The [scATAC-seq fragment files](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/output/fragments) are like bed files, but with a `.tsv` extension.

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

The fragments have been adjusted for the 9bp overhang created by Tn5 transposition so the ends are representative of the Tn5 cut sites. My analysis is primarily interested in inferring the exact Tn5 site that was bound, so we will convert the paired-end sequencing fragments to Tn5 sites.

This analysis uses 10 cell types that were used in the [ArchR](https://www.nature.com/articles/s41588-021-00790-6) paper. Each cluster of cells was identified and extracted as individual files containing the fragments associated with each cell type.

## Converting Fragments to Tn5 insertion sites

I convert the scATAC-seq fragment files to bed files, where each entry represents a Tn5 cut site that has been windowed to 40bp. I chose 40bp based on the x-ray crystallography of Tn5 that shows it has ~19bp wide binding site.

```
Fragment       ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓  
Tn5 IS         █                █
Slop 20bp    █████            █████
```

The 3 entries from the original fragments file are reported below:

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

