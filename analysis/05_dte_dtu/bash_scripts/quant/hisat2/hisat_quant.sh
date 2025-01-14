for bam in /Volumes/Fast_T7/htr8_test_output/test_bam/*.bam; do

    sample_name=$(basename ${bam} .bam)

    salmon quant --targets /Volumes/Fast_T7/Homo_sapiens.GRCh38.transcriptome.fa \
        --libType A \
        --alignments ${bam} \
        --output /Volumes/Fast_T7/htr8_test_output/test_quant/${sample_name}
done