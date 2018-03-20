#!/bin/bash

export TARGET=/n/desai_lab/users/ejerison/pleiotropy_sequencing_data/alignment_12_28_2017/sequencing_analysis_12_28_2017.tar.gz
export INPUT=/n/regal/desai_lab/ejerison/pleiotropy_sequencing/alignment_12_28_2017/sequencing_analysis_12_28_2017

sbatch --job-name=backup-$INPUT sequencing_analysis_12_28_2017/backup.sbatch