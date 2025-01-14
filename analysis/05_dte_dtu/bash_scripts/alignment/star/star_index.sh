#!/bin/bash

wd=/hpcfs/users/a1666761/290921_trophoblast_dtu
refs=/hpcfs/users/a1666761/Refs
genome=${refs}/ref_genomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa
annotation=${refs}/ref_annotations/Homo_sapiens.GRCh38.103.gtf

index=${refs}/STAR_index
cores=16

module load STAR/2.7.3a

# Generate STAR index
STAR --runThreadN ${cores} \
     --runMode genomeGenerate \
     --genomeDir ${index}/GRCh38 \
     --genomeFastaFiles ${genome} \
     --sjdbGTFfile ${annotation} \
     --sjdbOverhang 99
