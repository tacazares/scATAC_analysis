#!/usr/bin/python3

### Convert scATAC-seq fragments to a Tn5 bed file ###

# Inputs:

# input_fragments: Path to input fragment tsv file
# output directory: Path to directory to write results
# name: Output file name
# split: Split fragment ends into different files

# Outputs:

# A bed file of Tn5 sites that have been inferred from the scATAC-seq fragments file.
# If split is True: Then output two bed files, one for each end of the fragment.

######################################################

###-Imports-###
import pandas as pd
import argparse
import os

###-Arguments-###
parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_fragments",
                dest='FRAGMENTS', 
                help="Input fragments file",
                required=True)

parser.add_argument("-o", "--output_dir",
                dest='OUT_DIR', 
                help="Output directory",
                default="./output",
                required=False)

parser.add_argument("-n", "--name",
                dest='NAME', 
                help="Name to use for output file",
                default="Tn5_IS",
                required=False)

parser.add_argument("-split", "--split",
                action='store_true',
                dest='SPLIT', 
                help="Split fragment ends into seperate bed files",
                default="False",
                required=False)

args = parser.parse_args()

###-Functions-###
def convert_fragments_to_tn5_bed(fragments_tsv, chroms, split=False):
    """Convert scATAC fragments file to Tn5 insertion sites bed

    Args:
        fragments_tsv (str): Path to the scATAC-seq fragments file

    Returns:
        dataframe: Tn5 insertions sites in a BED formatted pandas dataframe
    """
    col_types = {
                 "chr": "category", 
                 "start": "int32", 
                 "stop": "int32", 
                 "barcode": "category", 
                 "support": "int32"
                }
    
    # Import fragments tsv as a dataframe
    df = pd.read_table(fragments_tsv,
                       header=None,
                       names=["chr", "start", "stop", "barcode", "support"],
                       dtype=col_types, 
                       low_memory=False)
    
    df = df[df["chr"].isin(chroms)]
    
    # subset starts
    df_starts = df[["chr", "start", "barcode"]].copy()
    
    # subset stops
    df_ends = df[["chr", "stop", "barcode"]].copy()
    
    # Add a 1 to create 1 bp intervals representing the cut site
    df_starts["stop"] = df_starts["start"].copy() + 1
    
    # Subtract a 1 bp interval to represent the cut site. Reference miralidlab wiki page for scATAC-seq analysis
    df_ends["start"] = df_ends["stop"].copy() - 1
        
    if split:
        # If split is True, return two dataframes, one for each end of the fragment    
        return df_starts[["chr", "start", "stop", "barcode"]], df_ends[["chr", "start", "stop", "barcode"]]
    
    else:
        # Concatenate both dataframes        
        df_cat = pd.concat([df_starts, df_ends])

        return df_cat[["chr", "start", "stop", "barcode"]]

###-Constants-###
ALL_CHRS = [
            "chr1", "chr2", "chr3", "chr4", "chr5", "chr6",
            "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
            "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
            "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"
            ]

###-Workflow-###
if args.SPLIT == True:
    # Create file names
    output_file_start = os.path.join(args.OUT_DIR, args.NAME + "_start.bed.gz")
    output_file_stop = os.path.join(args.OUT_DIR, args.NAME + "_stop.bed.gz")

    print("Splitting Fragment ends into individual files")
    
    start_interval_df, stop_interval_df = convert_fragments_to_tn5_bed(args.FRAGMENTS, ALL_CHRS, args.SPLIT)
    
    start_interval_df.to_csv(output_file_start, sep="\t", header=False, index=False)
    
    stop_interval_df.to_csv(output_file_stop, sep="\t", header=False, index=False)

else:
    # Create file names
    output_file = os.path.join(args.OUT_DIR, args.NAME + ".bed.gz")
    print("Splitting Fragment ends into a new BED file: " + output_file)

    tn5_bed = convert_fragments_to_tn5_bed(args.FRAGMENTS, ALL_CHRS)

    tn5_bed.to_csv(output_file, sep="\t", header=False, index=False)
