#!/bin/bash
#export population=$1
#export sequence_type=$2
export bam_files=""
export reference_file=~/w303_reference_genome/w303_ref.fasta
mkdir combined_mut_files

export sample_list=$(python sequencing_analysis_12_28_2017/file_name_utility.py sample_list)
for i in ${sample_list}; do
    export bam_files="${bam_files} breseq_output/bam_files/${i}.bam"
done
echo ${bam_files}

#rm data/timecourse_files/*
for i in $(python sequencing_analysis_12_28_2017/divide_chromosome.py keys); do
	export chrom=$i
	
	for j in $(python sequencing_analysis_12_28_2017/divide_chromosome.py $i)
#for i in {200000..250000..50000}
#for i in {50000..13000000..50000}
do
    export start_position=`expr $j - 50000`
    export end_position=$j
    #bash create_combined_muts_file_run.sh
    sbatch --job-name=merge_muts_${chrom}_${start_position} sequencing_analysis_12_28_2017/Step4_create_combined_muts_file.sbatch
	done
done
