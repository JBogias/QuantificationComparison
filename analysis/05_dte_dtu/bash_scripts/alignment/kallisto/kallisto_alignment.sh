#!/bin/bash -l

# Set paths to directories #
wd=/hpcfs/users/a1666761/290921_trophoblast_dtu
refs=/hpcfs/users/a1666761/Refs
output=${wd}/data/alignment/kallisto_global
data=${wd}/data/trimmed

# Run kallisto on list of samples
for fqgz in ${data}/*_R1*.fastq.gz; do

    SampleName=$(basename ${fqgz} _R1.fastq.gz)

    kallisto quant ${fqgz} ${fqgz/_R1/_R2} \
            -i ${refs}/kallisto_index/kallisto_transcripts_global.idx \
            -o ${output}/${SampleName} \
            -b 100 \
	    --threads 16
        
done
