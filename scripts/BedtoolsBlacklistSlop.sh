########################################################################
# Use bedtools to remove blacklisted regions and slop a BED file
########################################################################

###### INPUTS ######
bed_file=${1}
output_dir=${2}
blacklist_bed=${3}
chrom_sizes=${4}
slop_size=${5}

###### Parameters ######

# Set up output file names
slop_file=`basename ${bed_file} .bed.gz`_slop${slop_size}.bed.gz

echo "Making directory " ${output_dir}
mkdir -p ${output_dir}
cd ${output_dir}

###### WORKFLOW ######
echo "BED file: " ${bed_file} 
echo "Output directory: " ${output_dir} 
echo "Output File: " ${slop_size}
echo "Blacklist BED file: " ${blacklist_bed} 
echo "Chromosome sizes file: " ${chrom_sizes} 
echo "Slop size: " ${slop_size}

bedtools intersect -a ${bed_file} -b ${blacklist_bed} -v | bedtools slop -i - -g ${chrom_sizes} -b ${slop_size} | bedtools sort | pigz > ${slop_file}
