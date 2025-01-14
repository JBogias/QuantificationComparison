#!/bin/bash

ml arch/arch/skylake
ml SRA-Toolkit/2.10.8

dir=/hpcfs/users/a1666761/290921_trophoblast_dtu

# Because phoenix doesn't allow data downloading had to do it local
#prefetch --option-file ${dir}/files/SraAccList_dtu.txt

srrlist=${dir}/files/SraAccList_dtu.txt

while read srr;
  do
    fasterq-dump --outdir ${dir}/data/raw --split-files ${dir}/data/sra/${srr}/${srr}.sra
  done < ${srrlist}
