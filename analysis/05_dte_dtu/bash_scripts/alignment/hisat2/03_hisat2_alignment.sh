#!/bin/bash

wd=/path/to/wd
trim=${wd}/trimmed
idx=${wd}/index
sam=${wd}/sams


for fqgz in ${trim}/*_R1.fastq.gz; do
    SampleName=$(basename ${fqgz} _R1.fastq.gz)
    hisat2 -x ${idx}/grch38_genome_tran -1 ${fqgz} -2 ${fqgz/_R1/_R2} -S ${sam}/${SampleName}.sam
done