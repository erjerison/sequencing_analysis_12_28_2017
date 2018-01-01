#!/bin/bash
#SBATCH -p serial_requeue
#SBATCH -t 47:00:00
#SBATCH --mem 7400
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -o logs/create_part_muts.%j.out      #File to which standard out will be writte
#SBATCH -e logs/create_part_muts.%j.err      #File to which standard err will be written

export region="${chrom}:${start_position}-${end_position}"
echo $region
samtools mpileup -q10 -f ${reference_file} -r ${chrom}:${start_position}-${end_position} ${bam_files} > combined_mut_files/combined_muts_${chrom}_${start_position}.pileup

cat combined_mut_files/combined_muts_${chrom}_${start_position}.pileup | python sequencing_analysis_12_28_2017/create_combined_muts_file.py ${start_position} ${end_position} > combined_mut_files/combined_muts_${chrom}_${start_position}.txt

#bzcat realigned_timecourse.bz2 | python split_timecourse.py ${start_position} ${end_position} > data/timecourse_files/${start_position}_timecourse.txt 

#cat data/timecourse_files/${start_position}_timecourse.txt | ./annotate_pvalue > data/timecourse_files/${start_position}_pvalue_timecourse.txt  
