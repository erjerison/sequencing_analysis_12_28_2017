#!/bin/bash
#SBATCH -p serial_requeue
#SBATCH -t 47:00:00
#SBATCH --mem 32000
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -o logs/breseq.%j.out      #File to which standard out will be written
#SBATCH -e logs/breseq.%j.err      #File to which standard err will be written

mkdir -p breseq_output/${i}

#Call breseq once for each sample. In run 1 (plate 1) each sample was split across 2 lanes of a flow cell. For run 2 (plates 2, 3), each sample was in a single lane of the flow cell.
#For each lane, the reads for each sample have been split into the following files:
#1. unpaired reads that are the output of Trimmomatic. These are reads where only 1 read of the pair was retained because the reads completely overlapped.
#2. paired reads that are the output of Trimmomatic. Note that we are not merging overlapping paired-end reads because bases that disagree due to an error on one read or the other may be hard to call at this stage, but easier to call once we have piled up all reads that align to this base.

if [ -f trimmed_fastq_files/${i}*.gz ]; then
	for x in trimmed_fastq_files/${i}*.gz; do
		gunzip $x
	done
fi

breseq -p -j 1 -o breseq_output/${i} -r $REF_SEQ ${FASTQ_FILES} > breseq_output/${i}.out 2> breseq_output/${i}.err