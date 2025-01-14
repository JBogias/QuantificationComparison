#!/bin/bash

# Script to run RSEM quantification

# Packages needed (loaded in slurm script):
# arch/arch/haswell
# RSEM/1.2.30-foss-2017a

module load arch/arch/haswell
module load RSEM/1.2.30-foss-2017a

data=/hpcfs/users/a1666761/290921_trophoblast_dtu/data
refs=/hpcfs/users/a1666761/Refs
rsemref=${refs}/RSEM_index/GRCh38
cores=16

rsem-prepare-reference --gtf ${refs}/ref_annotations/Homo_sapiens.GRCh38.103.gtf \
                       ${refs}/ref_genomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
                       "${rsemref}/GRCh38_"
