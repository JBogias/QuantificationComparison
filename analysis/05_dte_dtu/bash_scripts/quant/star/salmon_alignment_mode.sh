#!/bin/bash

module load arch/arch/haswell
module load Salmon/1.1.0-foss-2016b

wd=/hpcfs/users/a1666761/290921_trophoblast_dtu
#transcriptome=/hpcfs/users/a1666761/Refs/ref_transcriptomes/gffread_grch38.fa
transcriptome=/hpcfs/users/a1666761/Refs/ref_transcriptomes/Homo_sapiens.GRCh38.transcriptome.fa

for bam in ${wd}/data/alignment/star/*toTranscriptome.out.bam; do

    sample_name=$(basename ${bam} _GRCh38_Aligned.toTranscriptome.out.bam)

    salmon quant --targets ${transcriptome} \
        --libType ISR \
	--numBootstraps 100 \
        --alignments ${bam} \
        --output ${wd}/data/quant/star_salmon_quant/${sample_name}

done
