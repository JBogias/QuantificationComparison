#!/bin/bash

wd=/path/to/wd
sam=${wd}/sams
bam=${wd}/bams

for sambam in ${sam}/*.sam; do
    SampleName=$(basename ${sambam} .sam)
    samtools view -bS ${sam} > ${bam}/${SampleName}.bam
done