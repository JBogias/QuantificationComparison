#!/bin/bash

# Script to run STAR

# Packages needed:
#	STAR

dir=/hpcfs/users/a1666761/290921_trophoblast_dtu
refs=/hpcfs/users/a1666761/Refs
data=${dir}/data/trimmed
output=${dir}/data/alignment/star

index=${refs}/STAR_index/GRCh38

build=GRCh38
cores=16

module load STAR

# Align trimmed data to reference genome
for FQGZ in ${data}/*_R1*.fastq.gz; do

	SampleName=$(basename ${FQGZ} _R1.fastq.gz)

	STAR --genomeDir ${index} \
		--readFilesIn ${FQGZ} ${FQGZ/_R1/_R2} \
		--readFilesCommand zcat \
		--outFilterType BySJout \
		--alignSJoverhangMin 7 \
		--outFilterMismatchNmax 999 \
		--outSAMtype BAM Unsorted \
		--quantMode TranscriptomeSAM \
		--outFileNamePrefix ${output}/"${SampleName}"_"${build}"_ \
		--outSAMattrRGline ID:"${SampleName}" LB:library PL:illumina \
				   PU:machine SM:"${build}" \
		--outSAMmapqUnique 60 \
		--runThreadN ${cores}
done
