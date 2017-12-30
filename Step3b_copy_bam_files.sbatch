#!/bin/bash
#SBATCH -p serial_requeue
#SBATCH -t 2:00:00
#SBATCH --mem 3750
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -o logs/copy.%j.out      #File to which standard out will be written
#SBATCH -e logs/copy.%j.err      #File to which standard err will be written

cp $old_path_bam $new_path_bam