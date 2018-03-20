#!/bin/bash
export Array=(chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chr12 chr13 chr14 chr15 chr16)
for i in ${Array[*]}; do
	export chrom=$i
	export new_file=mutation_data/all_muts_alts_covs.txt
    export filename=mutation_data/${chrom}_potential_de_novo_muts.txt
	tail -n +2 ${filename} >> ${new_file}
	
done
