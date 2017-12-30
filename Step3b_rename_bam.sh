#!bin/bash/
mkdir breseq_output/bam_files
#This is a utility for putting all of the alignment files output by Breseq into one folder
for i in $(python sequencing_analysis_12_28_2017/file_name_utility.py sample_list); do
	old_path_index=breseq_output/${i}/data/reference.bam.bai
	new_path_index=breseq_output/bam_files/${i}.bam.bai
	old_path_bam=breseq_output/${i}/data/reference.bam
	new_path_bam=breseq_output/bam_files/${i}.bam
	echo "$old_path_bam"
	export old_path_bam
	export new_path_bam
	cp $old_path_index $new_path_index
	sbatch --job-name=copy_bam_files-${i} sequencing_analysis_12_28_2017/Step3b_copy_bam_files.sbatch
done

