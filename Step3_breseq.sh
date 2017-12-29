#!/bin/bash

mkdir -p breseq_output

#for i in $(python  sequencing_analysis_12_28_2017/file_name_utility.py sample_list); do
#for i in ${"E1-1"}; do
    export i=E1-1
    echo $i
    export FASTQ_FILES=$(python sequencing_analysis_12_28_2017/file_name_utility.py $i)
    export REF_SEQ=~/S288c_ref.gff
    #export REF_SEQ=~/w303_reference_genome/w303_ref_with_seq_no_2micron.gff3
    sbatch --job-name=breseq-${i} sequencing_analysis_12_28_2017/Step3_breseq.sbatch
#done