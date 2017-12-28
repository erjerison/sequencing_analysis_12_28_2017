#!/bin/bash
mkdir trimmed_fastq_files

for i in $(python sequencing_analysis_12_28_2017/file_name_utility.py trial); do
    export i
    echo $i
    sbatch --job-name=trim-${i} sequencing_analysis_12_28_2017/Step2_trim_adapters.sbatch
done