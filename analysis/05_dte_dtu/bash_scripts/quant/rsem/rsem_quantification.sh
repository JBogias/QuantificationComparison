#!/bin/bash

# Script to run RSEM quantification

# Packages needed (loaded in slurm script):
# arch/arch/haswell
# RSEM/1.2.30-foss-2017a

data=/hpcfs/users/a1666761/210316_GDM_RNAseq/data/alignment/STAR_txp
refs=/hpcfs/users/a1666761/Refs
rsemref=${refs}/RSEM_index
output=/hpcfs/users/a1666761/210316_GDM_RNAseq/data/quant/RSEM
cores=16

rsem-prepare-reference --gtf ${refs}/ref_annotations/Homo_sapiens.GRCh37.87.gtf \
                       ${refs}/ref_genomes/Homo_sapiens.GRCh37.dna.primary_assembly.fa \
                       ${rsemref}/GRCh37_genome

for BAM in ${data}/*_GRCh37_Aligned.toTranscriptome.out.bam; do

    SampleName=$(basename ${BAM} _GRCh37_Aligned.toTranscriptome.out.bam)

    rsem-calculate-expression --paired-end \
			      --no-bam-output \
			      --quiet \
			      --seed 12345 \
			      --num-threads ${cores} \
	  		      --alignments \
			      ${BAM} \
                              ${rsemref}/GRCh37_genome \
			      ${output}/${SampleName}

done
