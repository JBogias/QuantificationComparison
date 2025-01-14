#!/bin/bash

wd=/path/to/wd
idx=${wd}/index
fqs=${wd}/fastq
trim=${wd}/trimmed
sam=${wd}/sams
bam=${wd}/bams
refs=${wd}/refs
quant=${wd}/quant

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

hisat2-build -p 16 --exon ${refs}/hs_grch38_exons.txt \
    --ss ${refs}/hs_grch38_splice_sites.txt \
    ${refs}/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    ${idx}/grch38_genome_tran

for fqgz in ${trim}/*_R1.fastq.gz; do
    SampleName=$(basename ${fqgz} _R1.fastq.gz)
    hisat2 -x ${idx}/grch38_genome_tran -1 ${fqgz} -2 ${fqgz/_R1/_R2} -S ${sam}/${SampleName}.sam
done

for sambam in ${sam}/*.sam; do
    SampleName=$(basename ${sambam} .sam)
    samtools view -bS ${sam} > ${bam}/${SampleName}.bam
done

for bamfile in ${bam}/*.bam; do
    SampleName=$(basename ${bam} .bam)
    salmon quant --targets ${refs}/Homo_sapiens.GRCh38.transcriptome.fa \
        --libType A \
        --alignments ${bamfile} \
        --output ${quant}/${sample_name}
done

