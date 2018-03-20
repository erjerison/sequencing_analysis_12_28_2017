#!/bin/bash
#SBATCH -p serial_requeue
#SBATCH -t 47:00:00
#SBATCH --mem 37500
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -o logs/backup.%j.out      #File to which standard out will be writte
#SBATCH -e logs/backup.%j.err      #File to which standard err will be written

tar -cvzf $TARGET $INPUT