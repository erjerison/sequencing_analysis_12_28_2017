import numpy
import sys
from file_name_utility2 import sample_list

samples = sample_list()

input_file = 'mutation_data/all_muts_alts_covs_annotated_with_aa_positions.txt'
output_file = 'mutation_data/mutation_lists_with_aa_positions_1_3_2018.txt'

file = open(input_file,'r')

mut_dict = {}
gene_name_column = 2*253 + 4

hap_dip_file = open('/Users/ejerison/Dropbox/Pleiotropy_Experiment/paper_analyses/data/ploidy_screen_results.txt', 'r')

hap_dip_dict = {}

firstline = True
for line in hap_dip_file:
	if firstline:
		firstline = False
	else:
	
		line_list = line.strip().split(', ')
		env,clone = line_list[0].split('_')
		new_label = 'E' + env.strip('evolenv') + '-' + clone.strip('clone')
		
		hap_dip_dict[ new_label ] = line_list[1]

hap_dip_file.close()
print(hap_dip_dict)
first_line = 1
for line in file:
	
	alts_covs = line.strip().split(',')
	alts = [float(a) for a in alts_covs[4].split()[0:253]]
	covs = [float(a) for a in alts_covs[5].split()[0:253]]
	
	mutation = alts_covs[0:4]
	mutation.extend(alts_covs[5].split()[253:])
	mut_presence = numpy.array(numpy.array(alts) > 1.5)
	
	sample_idx = numpy.where(mut_presence > .5)[0]
	
	for entry in sample_idx:
		freq = alts[entry]/covs[entry]
		
		if samples[entry] in hap_dip_dict and hap_dip_dict[ samples[entry] ] == 'D': ##This is a diploid, so we will call hets as well as homozygotes
		
			if (samples[entry] in mut_dict and freq > .4 and freq < .6 and covs[entry] > 10):
				mut_dict[samples[entry]].append((mutation, 'Het'))
			elif (freq > .4 and freq < .6 and covs[entry] > 10):
				mut_dict[samples[entry]] = [(mutation, 'Het')]
			elif (samples[entry] in mut_dict and freq > .8):
				mut_dict[samples[entry]].append((mutation, 'Hom'))
			elif freq > .8:
				mut_dict[samples[entry]] = [(mutation, 'Hom')]
		else:
			if (samples[entry] in mut_dict and freq > .8):
				mut_dict[samples[entry]].append((mutation, 'NA'))
			elif freq > .8:
				mut_dict[samples[entry]] = [(mutation, 'NA')]
		
file.close()


		
#print mut_dict
	
file2 = open(output_file,'w')
for sample in sorted(samples):
	
	file2.write(sample + '\n')
	if sample in mut_dict:
		for entry in mut_dict[sample]:
			mutation_strs = ('\t').join(str(e) for e in entry[0])
			entry_strs = (',').join((mutation_strs,str(entry[1])))
			file2.write(entry_strs + '\n')
file2.close()
		
