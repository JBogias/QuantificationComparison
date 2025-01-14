#!/bin/bash

# Trimming of HTR8/svneo cell data from Georgiadou et al. 2021

dir=/hpcfs/users/a1666761/290921_trophoblast_dtu
data=${dir}/data/raw
output=${dir}/data/trimmed

for fastqgz in ${data}/*sra_1.fastq.gz; do

    sample_name=$(basename ${fastqgz} .sra_1.fastq.gz)

    AdapterRemoval --file1 ${fastqgz} \
                   --file2 ${fastqgz/sra_1/sra_2} \
                   --output1 ${output}/${sample_name}_R1.fastq.gz \
                   --output2 ${output}/${sample_name}_R2.fastq.gz \
                   --trimns \
                   --trimqualities \
                   --threads 16 \
                   --gzip
done

