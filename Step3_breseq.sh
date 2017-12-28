#!/bin/bash

mkdir -p breseq_output

#for i in $(python scripts/file_name_utility.py sample_list); do
for i in ${E1-1}; do
    export i
    echo $i
    export FASTQ_FILES=$(python sequencing_analysis_12_28_2017/file_name_utility.py $i)
    export REF_SEQ=~/w303_reference_genome/w303_ref_with_seq.gff
    sbatch --job-name=breseq-${i} scripts/Step3_breseq.sbatch
done
