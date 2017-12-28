#!/bin/bash
#SBATCH -p serial_requeue
#SBATCH -t 5:00:00
#SBATCH --mem-per-cpu 16000
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -o logs/trim.out      #File to which standard out will be written
#SBATCH -e logs/trim.err      #File to which standard err will be written

#Because of the bzipped files, we will unzip before trimming
##Some of the files are gzipped
for x in fastqs_with_sample_names/${i}*.gz; do
	if [ -e "$x" ]; then
		gunzip $x
	fi
done

##Some of the files are bzipped
for x in fastqs_with_sample_names/${i}*.bz2; do
	if [ -e "$x" ]; then
		bzip2 -d $x
	fi
done

# adaptor trimming
java -jar $TRIMMOMATIC/trimmomatic-0.32.jar PE -threads 1 -phred33 fastqs_with_sample_names/${i}.R1.fastq fastqs_with_sample_names/${i}.R2.fastq trimmed_fastq_files/${i}.R1.paired.gz trimmed_fastq_files/${i}.R1.unpaired.gz trimmed_fastq_files/${i}.R2.paired.gz trimmed_fastq_files/${i}.R2.unpaired.gz ILLUMINACLIP:$TRIMMOMATIC/adapters/NexteraPE-PE.fa:2:30:10:2:'false' LEADING:20 TRAILING:20 2> trimmed_fastq_files/${i}.trim.txt
#java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -threads 1 -phred33 fastqs_with_sample_names/${i}.R1.fastq fastqs_with_sample_names/${i}.R2.fastq trimmed_fastq_files/${i}.R1.paired.gz trimmed_fastq_files/${i}.R1.unpaired.gz trimmed_fastq_files/${i}.R2.paired.gz trimmed_fastq_files/${i}.R2.unpaired.gz ILLUMINACLIP:Trimmomatic-0.35/adapters/NexteraPE-PE.fa:2:30:10:2:'false' LEADING:20 TRAILING:20 2> trimmed_fastq_files/${i}.trim.txt