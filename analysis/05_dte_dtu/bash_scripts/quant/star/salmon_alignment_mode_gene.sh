#!/bin/bash

module load arch/arch/haswell
module load Salmon/1.1.0-foss-2016b

wd=/hpcfs/users/a1666761/290921_trophoblast_dtu
transcriptome=/hpcfs/users/a1666761/Refs/ref_genomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa
#transcriptome=/hpcfs/users/a1666761/Refs/ref_transcriptomes/Homo_sapiens.GRCh38.transcriptome.fa

for bam in ${wd}/data/alignment/star/*_sortedByCoordinate.out.bam; do

    sample_name=$(basename ${bam} __sortedByCoordinate.out.bam)

    salmon quant --targets ${transcriptome} \
        --libType ISR \
        --alignments ${bam} \
        --output ${wd}/data/quant/star_gene_level_salmon/${sample_name}

done
