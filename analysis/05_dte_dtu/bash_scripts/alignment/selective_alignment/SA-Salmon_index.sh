#!/bin/bash -l

# Script by Justin Bogias #

# Purpose: Perform selective alignment followed by Salmon quantification on a
#  list of fastq files. This script relies on decoy transcripts which are
#  generated using a script obtained from the SalmonTools github repository

# Set variables

output=/hpcfs/users/a1666761/290921_trophoblast_dtu/data
transcripts=/hpcfs/users/a1666761/Refs/ref_transcriptomes/GRCh38_decoy_global/gentrome.fa
decoys=/hpcfs/users/a1666761/Refs/ref_transcriptomes/GRCh38_decoy_global/decoys.txt

# Required packages
module load arch/arch/haswell
module load Salmon/1.1.0-foss-2016b

# Index with decoy transcriptome
salmon index -t ${transcripts} \
             -i ${output}/indexes/salmon_decoy_index_GRCh38 \
             --decoys ${decoys} \
             -k 31
