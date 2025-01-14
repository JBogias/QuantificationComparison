#!/bin/bash -l

# Script to run samtools sort

# Packages needed:
#       SAMtools

data=/hpcfs/users/a1666761/290921_trophoblast_dtu/data/alignment/star
cores=16

for BAM in ${data}/*_Aligned_sortedByCoordinate.out.bam; do

    SampleName=$(basename ${BAM} .bam)

    samtools index ${BAM} ${data}/${SampleName}.bai

done

