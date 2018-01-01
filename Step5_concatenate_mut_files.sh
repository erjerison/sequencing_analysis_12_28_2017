#!/bin/bash

mkdir combined_mut_files/full_chromosomes
for i in $(python sequencing_analysis_12_28_2017/divide_chromosome.py keys); do
	export chrom=$i
	export new_file=combined_mut_files/full_chromosomes/combined_muts_${chrom}.txt
	rm ${new_file}
	
	for j in $(python sequencing_analysis_12_28_2017/divide_chromosome.py $i)

do
    export start_position=`expr $j - 50000`
    export filename=combined_mut_files/combined_muts_${chrom}_${start_position}.txt
	tail -n +2 ${filename} >> ${new_file}
	done
done
