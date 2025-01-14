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

transcripts=${refs}/ref_transcriptomes/GRCh38_decoy_global/gentrome.fa
decoys=${refs}/ref_transcriptomes/GRCh38_decoy_global/decoys.txt

# Index with decoy transcriptome

# salmon index -t ${transcripts} \
#              -i ${output}/indexes/salmon_decoy_index_GRCh38 \
#              --decoys ${decoys} \
#              -k 31

index=/hpcfs/users/a1666761/290921_trophoblast_dtu/data/indexes/salmon_decoy_index_GRCh38

# Run Selective alignment and Salmon quant on list of samples
for sample in ${data}/*R1.fastq.gz;
do

    sample_name=$(basename ${sample} _R1.fastq.gz)

    salmon quant -i ${index} \
           --libType A \
           -1 ${sample} \
           -2 ${sample/_R1/_R2} \
           --numBootstraps 0 \
           --validateMappings \
           --dumpEq \
           --threads 16 \
           -o ${output}/${sample_name}_b01

done
