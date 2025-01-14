#!/bin/bash

wd=/path/to/wd
bam=${wd}/bams
refs=${wd}/refs
quant=${wd}/quant

for bamfile in ${bam}/*.bam; do
    SampleName=$(basename ${bam} .bam)
    salmon quant --targets ${refs}/Homo_sapiens.GRCh38.transcriptome.fa \
        --libType A \
        --alignments ${bamfile} \
        --output ${quant}/${sample_name}
done