#!/bin/bash -l

# Script to run samtools sort

# Packages needed:
#	SAMtools

data=/hpcfs/users/a1666761/290921_trophoblast_dtu/data/alignment/star
cores=16

for bam in ${data}/*_Aligned.out.bam; do

	SampleName=$(basename ${bam} Aligned.out.bam)

	samtools sort ${bam} \
		 -o ${data}/${SampleName}_sortedByCoordinate.out.bam \
		 --threads=${cores}

done
