#!/bin/bash
mkdir trimmed_fastq_files

for i in $(python scripts/file_name_utility.py trial); do
    export i
    echo $i
    sbatch --job-name=trim-${i} scripts/Step2_trim_adapters.sbatch
done
