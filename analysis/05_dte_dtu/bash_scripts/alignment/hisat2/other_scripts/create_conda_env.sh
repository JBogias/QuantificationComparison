conda create -n hisat2_rnaseq

conda install -c bioconda salmon=1.10.1
conda install -c bioconda samtools=1.6
conda install -c bioconda AdapterRemoval=2.2.2
conda install -c bioconda hisat2=2.2.1

conda deactivate