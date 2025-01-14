#!/bin/bash -l

# Script by Justin Bogias #

# Purpose: Perform selective alignment followed by Salmon quantification on a
#  list of fastq files. This script relies on decoy transcripts which are
#  generated using a script obtained from the SalmonTools github repository

# Set variables

# Chromosome patches, supercontigs and haplotypes were removed from the assembly
# grep -vE  "(HG|HSC|GL)" gentrome.fa > gentrome_nopatch.fa
# grep -vE  "(HG|HSC|GL)" decoys.txt > decoys_nopatch.txt

path=/hpcfs/users/a1666761/290921_trophoblast_dtu
output=${path}/data/alignment/selective_alignment
data=${path}/data/trimmed

info=${path}/files/SraAccList_dtu.txt
refs=/hpcfs/users/a1666761/Refs

transcripts=${refs}/ref_transcriptomes/GRCh38_decoy/gentrome.fa
decoys=${refs}/ref_transcriptomes/GRCh38_decoy/decoys.txt

# Index with decoy transcriptome

# salmon index -t ${transcripts} \
#              -i ${output}/indexes/salmon_decoy_index_GRCh38 \
#              --decoys ${decoys} \
#              -k 31

index=${output}/indexes/salmon_decoy_index_GRCh38

# Run Selective alignment and Salmon quant on list of samples
sample1=/hpcfs/users/a1666761/290921_trophoblast_dtu/data/trimmed/SRR13401123_R1.fastq.gz
sample2=/hpcfs/users/a1666761/290921_trophoblast_dtu/data/trimmed/SRR13401123_R2.fastq.gz

    sample_name=$(basename ${sample1} _R1.fastq.gz)

    salmon quant -i ${index} \
           --libType A \
           -1 ${sample1} \
           -2 ${sample2} \
           --numBootstraps 100 \
           --validateMappings \
           --dumpEq \
           --threads 16 \
           -o ${output}/${sample_name}
