#!/bin/bash

wd=/hpcfs/users/a1666761/290921_trophoblast_dtu
refs=/hpcfs/users/a1666761/refs
data=${wd}/data/trimmed
output=${wd}/data/alignment/star
counts=${wd}/data/featureCounts

genome=${refs}/ref_genomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
annotation=${refs}/ref_annotation/Homo_sapiens.GRCh38.103.gtf.gz

index=${refs}/indexes

module load STAR

build="GRCh38"
cores=16

# Generate STAR index
STAR --runThreadN ${cores} \
     --runMode genomeGenerate \
     --genomeDir ${index}/STAR_genome_index \
     --genomeFastaFiles ${genome} \
     --sjdbGTFfile ${annotation} \
     --sjdbOverhang 99 \
     --limitGenomeGenerateRAM 90043555797

# Align trimmed data to reference genome
for fqgz in ${data}/*R1*.fastq.gz; do

	SampleName=$(basename ${fqgz} _R1.fastq.gz)

	STAR --genomeDir ${index} \
	     --readFilesIn ${fqgz} ${fqgz/R1/R2} \
	     --readFilesCommand zcat \
	     --outFilterType BySJout \
	     --outFilterMismatchNmax 999 \
	     --outSAMtype BAM Sorted \
	     --outFileNamePrefix ${output}/"${SampleName}"_"${build}"_ \
	     --outSAMattrRGline ID:"${SampleName}" \
				LB:library \
				PL:illumina \
				PU:machine \
				SM:"${build}" \
	     --outSAMmapqUnique 60 \
	     --runThreadN ${cores}
done

# Now run featuresCounts
for bams in ${output}/*.bam.gz; do
	featureCounts -a ${annotation} \
		      -g "gene_id" \
		      -Q 10 \
		      -o ${counts}/gene_counts_quality.tsv \
		      -T 16 \
		      -s 1 ${output}/${bams}.bam.gz
done
