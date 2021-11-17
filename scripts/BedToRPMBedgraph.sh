########################################################################
# Convert a BED file into a library normalized bedgraph
########################################################################

# This code is for those interested in generating only the bedgraph, but not compressing it. Mainly for testing purposes

###### INPUTS ######
bed_file=${1}
output_dir=${2}
chrom_sizes=${3}

echo ${bed_file} ${output_dir} ${chrom_sizes}

###### Parameters ######

# Set up output file names
bedgraph_file=`basename ${bed_file} .bed.gz`.bg

echo "Making directory " ${output_dir}
mkdir -p ${output_dir}
cd ${output_dir}

###### WORKFLOW ######
echo "Calculate scaling factors for genome coverage"
read_count=$(gunzip -c ${bed_file} | wc -l)

# Generate bedgraph of cutsite coverage and normalize to sequencing depth
reads_factor=$(bc -l <<< "1/${read_count}")
rpm_factor=$(bc -l <<< "${reads_factor} * 1000000")

echo "Scaling factor " ${rpm_factor}

echo "Generate coverage track that is RPM normalized"
bedtools genomecov -i ${bed_file} -g ${chrom_sizes} -bg -scale ${rpm_factor} > ${output_dir}/${bedgraph_file}

echo "Done!"