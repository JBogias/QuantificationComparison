#!/bin/bash -l

# Set paths to directories #
wd=/hpcfs/users/a1666761/290921_trophoblast_dtu
refs=/hpcfs/users/a1666761/Refs
output=${wd}/data/alignment/kallisto
data=${wd}/data/trimming

transcripts=${refs}/ref_transcriptomes/Homo_sapiens.GRCh38.transcriptome.fa

# Generate index
kallisto index -i ${refs}/kallisto_index/kallisto_transcripts_global.idx ${transcripts}

