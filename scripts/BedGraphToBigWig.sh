#!/bin/bash

########################################################################
# Convert a BedGraph file into a BigWig
########################################################################

input_bedgraph=${1}
output_bigwig=${2}
chrom_sizes=${3}

bedGraphToBigWig ${1} ${3} ${2}
