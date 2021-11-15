########################################################################
# Convert a Tn5 cut site file into a library normalized bigwig
########################################################################

###### INPUTS ######
slop_file=${1}
output_dir=${2}
blacklist_bed=${3}
chrom_sizes=${4}

echo ${slop_file} ${output_dir} ${blacklist_bed} ${chrom_sizes}

###### Parameters ######

# Set up output file names
bedgraph_file=`basename ${slop_file} .bed`.bg
bigwig_file=`basename ${slop_file} .bed`.bw

echo "Making directory " ${output_dir}
mkdir -p ${output_dir}
cd ${output_dir}

###### WORKFLOW ######
echo "Calculate scaling factors for genome coverage"
read_count=$(wc -l < ${slop_file})

# Generate bedgraph of cutsite coverage and normalize to sequencing depth
reads_factor=$(bc -l <<< "1/${read_count}")
rpm_factor=$(bc -l <<< "${reads_factor} * 1000000")

echo "Scaling factor " ${rpm_factor}

echo "Generate coverage track that is RPM normalized"
bedtools intersect -a ${slop_file} -b ${blacklist_bed} -v | \
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