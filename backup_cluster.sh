#!/bin/bash

export TARGET=blah
export INPUT=bloo

sbatch --job-name=backup-$INPUT sequencing_analysis_12_28_2017/backup.sbatch