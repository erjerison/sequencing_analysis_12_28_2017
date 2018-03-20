import numpy
import sys
from file_name_utility2 import sample_list
import matplotlib.pylab as pt

samples = sample_list()

input_file = 'mutation_data/all_muts_alts_covs_annotated_with_aa_positions.txt'
output_file = 'mutation_data/coverage_by_clone.txt'
header_file = 'mutation_data/chr01_potential_de_novo_muts.txt'

file = open(header_file, 'r')
line = file.readline()
file.close()

clone_list = line.split(',')[4].split()

file = open(input_file,'r')

mut_dict = {}
gene_name_column = 2*253 + 4

n_lines = 0
net_covs = numpy.zeros((253,))

first_line = 1
for line in file:
	
	alts_covs = line.strip().split(',')
	alts = [float(a) for a in alts_covs[4].split()[0:253]]
	covs = [float(a) for a in alts_covs[5].split()[0:253]]
	
	net_covs += numpy.array(covs)
	n_lines += 1
	
file.close()
avg_cov = net_covs/float(n_lines)

print("average coverage across all clones,", numpy.mean(avg_cov))

file_out = open(output_file, 'w')
file_out.write( ("\t").join(clone_list) + "\n")
file_out.write( ("\t").join([str(cov) for cov in avg_cov]))
file_out.close()

#pt.hist(avg_cov, bins = numpy.arange(40))

for i in range(253):
	if avg_cov[i] < 5:
		print(clone_list[i], avg_cov[i])

#pt.show()
