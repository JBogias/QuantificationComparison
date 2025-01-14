#!/bin/bash

wd=/hpcfs/users/a1666761/290921_trophoblast_dtu
data=${wd}/data/trimmed
index=/hpcfs/users/a1666761/Refs/STAR_index
output=${wd}/data/alignment/star

build="GRCh38"
cores=16

module load STAR/2.7.3a

# Align trimmed data to reference genome
for fqgz in ${data}/*_R1*.fastq.gz; do

        SampleName=$(basename ${fqgz} _R1.fastq.gz)

        STAR --genomeDir ${index}/GRCh38 \
             --readFilesIn ${fqgz} ${fqgz/_R1/_R2} \
             --readFilesCommand zcat \
             --outFilterType BySJout \
             --outFilterMismatchNmax 999 \
             --outSAMtype BAM Unsorted \
             --outFileNamePrefix ${output}/"${SampleName}"_"${build}"_ \
             --outSAMattrRGline ID:"${SampleName}" \
                                LB:library \
                                PL:illumina \
                                PU:machine \
                                SM:"${build}" \
             --outSAMmapqUnique 60 \
             --runThreadN ${cores}
done
