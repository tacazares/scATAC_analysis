  #!/bin/bash

############### sra2fastq_scATAC.sh ###############
# This script will download and compress a fastq file from the SRA database

# Inputs:
# ${1}: SRR ID

# Outputs:
# gzipped fastq file

# Workflow:
# 1) download
# 2) convert to fastq and compress with gzip
##################################################

# Use wget to download the SRA from the link
prefetch ${1}

# Use fasterq-dump to split a scRNA SRA file into three fastq files
##fasterq-dump -e 8 --split-files ${SRA} 
fastq-dump --origfmt --gzip --split-files ${1}
