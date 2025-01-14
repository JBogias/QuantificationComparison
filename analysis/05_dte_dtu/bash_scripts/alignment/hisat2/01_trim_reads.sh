#!/bin/bash

wd=/path/to/wd
fqs=${wd}/fastq
trim=${wd}/trimmed


for fastqgz in ${fqs}/*sra_1.fastq.gz; do
    SampleName=$(basename ${fastqgz} .sra_1.fastq.gz)
    AdapterRemoval --file1 ${fastqgz} \
                   --file2 ${fastqgz/sra_1/sra_2} \
                   --output1 ${trim}/${SampleName}_R1.fastq.gz \
                   --output2 ${trim}/${SampleName}_R2.fastq.gz \
                   --trimns \
                   --trimqualities \
                   --threads 16 \
                   --gzip
done