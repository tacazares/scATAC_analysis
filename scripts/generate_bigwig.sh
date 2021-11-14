########################################################################
# Convert a Tn5 cut site file into a library normalized bigwig
########################################################################

###### INPUTS ######
output_dir=${2}
mapped_reads=${3}
blacklist_bed=${4}
chrom_sizes=${5}

insertion_sites
bedgraph_file=`basename ${slop_file} .bed`.bg
bigwig_file=`basename ${bedgraph_file} .bg`.bw

###### WORKFLOW ######
echo "Making directory " ${output_dir}
mkdir -p ${output_dir}
cd ${output_dir}

echo "Calculate scaling factors for genome coverage with bedtools"
# Generate bedgraph of cutsite coverage and normalize to sequencing depth
reads_factor=$(bc -l <<< "1/${mapped_reads}")
rpm_factor=$(bc -l <<< "${reads_factor} * 1000000")
bedgraph_file=`basename ${slop_file} .bed`.bg
out_filename_bw=`basename ${slop_file} .bed`.bw

echo "Scaling factor " ${rpm_factor}

echo "Generate coverage track that is RPM normalized"
bedtools intersect -a ${output_dir}/${slop_file} -b ${blacklist_bed} -v | \
bedtools sort | \
bedtools genomecov -i - -g ${chrom_sizes} -bg -scale ${rpm_factor} > \
${output_dir}/${bedgraph_file}

echo "Use bedGraphToBigWig to generate bw file"
bedGraphToBigWig ${output_dir}/${bedgraph_file} \
${chrom_sizes} \
${output_dir}/${bigwig_file}

echo "Cleanup"
rm ${output_dir}/${bedgraph_file}

echo "Done!"