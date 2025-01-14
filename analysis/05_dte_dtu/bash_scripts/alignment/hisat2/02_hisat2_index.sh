#!/bin/bash

wd=/path/to/wd
idx=${wd}/index
refs=${wd}/refs

# Make the splice and exon file using HISAT2 scripts
# extract_splice_sites.py ${refs}/Homo_sapiens.GRCh38.103.gtf > ${refs}/hs_grch38_splice_sites.txt
# extract_exons.py ${refs}/Homo_sapiens.GRCh38.103.gtf > ${refs}/hs_grch38_exons.txt

hisat2-build -p 16 --exon ${refs}/hs_grch38_exons.txt \
    --ss ${refs}/hs_grch38_splice_sites.txt \
    ${refs}/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    ${idx}/grch38_genome_tran