import numpy
import sys
import matplotlib.pylab as pt

import random
import string
complem = string.maketrans('GATC','CTAG')

##First read in the list of called mutations to check for cases where there are two snps, or a snp and an indel, in the same codon.
##These will be mis-annotated.

##Read in the genome and locations of genes
## read in the genetic code dictionary
file = open('genetic_code_table.txt','r')

genetic_code = {}
for line in file:
	linelist = line.strip().split('\t')
	codons = linelist[2].split(', ')
	for i in range(len(codons)):
		genetic_code[codons[i]] = linelist[0].strip()
file.close()		

numsamples = 253

chrom_labels = [''.join(('chr',i)) for i in ['I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII','XIV','XV','XVI']]
chrom_dict = dict(zip(chrom_labels,range(16)))

sequence = {}

read_list = []

gene_list = []
gene_start_stop = {}
gene_length = {}
strand = {}
phase = {}
gene_list_file = 'w303_ref_with_seq.gff'

file_genes = open(gene_list_file,'r')

file_line_list = list(file_genes)

for i in range(len(file_line_list)):
	if file_line_list[i].startswith('#'):
		continue
	elif file_line_list[i].startswith('>'):
		break
	else:
		linelist_prev = file_line_list[i-1].split('\t')
		linelist = file_line_list[i].split('\t')
		#print linelist
		if linelist[2] == 'gene':
			#record the gene name, whether it is
			gene_name = linelist[8].strip().split(';')[1].split('=')[1]
			entry_label = linelist[8].strip().split(';')[1].split('=')[0] 
			
		if linelist[2] == 'CDS':
			#record the location and phase of particular CDS's within the gene
			
			parent = linelist[8].strip().split(';')[1].split('=')[1].split('_')[0]
			
			if parent == gene_name:
				if entry_label != 'gene':
					gene_list.append(parent)
					gene_start_stop[parent] = [float(linelist[3]),float(linelist[4])]
					gene_length[parent] = float(linelist[4]) - float(linelist[3])
					strand[parent] = linelist[6]
					phase[parent] = linelist[7]
				else:
					gene_list.append(gene_name)
					gene_start_stop[gene_name] = [float(linelist[3]),float(linelist[4])]
					gene_length[gene_name] = float(linelist[4]) - float(linelist[3])
					strand[gene_name] = linelist[6]
					phase[gene_name] = linelist[7]
			else:
				name = linelist[8].strip().split(';')[2].split('=')[1]
				gene_list.append(name)
				gene_start_stop[name] = [float(linelist[3]),float(linelist[4])]
				gene_length[name] = float(linelist[4]) - float(linelist[3])
				
				strand[name] = linelist[6]
				phase[name] = linelist[7]

num_features = len(gene_start_stop)

genome_location = 0

for j in numpy.arange(i,len(file_line_list)):
	if file_line_list[j].startswith('>'):
		chr = file_line_list[j].strip()[1:]
		sequence[chr] = []
	
	elif len(sequence[chr]) == 0:
		sequence[chr].append(file_line_list[j].strip())
	else:
		sequence[chr][0] += file_line_list[j].strip()

file_genes.close()

input_file = 'mutation_data/mutation_lists_with_aa_positions_1_3_2018.txt'

#First loop to record the amino acid of all mutations in genes

file = open(input_file,'r')
file_lines = file.readlines()
file.close()

amino_acid_dict = {}
mutation_dict = {}
for line in file_lines:
	linelist = line.strip().split('\t')
	if len(linelist) < 1.5:
		#Go to the next clone
		clone_name = linelist[0]
		amino_acid_dict[clone_name] = {}
		mutation_dict[clone_name] = {}
	else:
		mutation = ('_').join(str(i) for i in linelist)
		if len(linelist) > 5.5:
			gene_name = linelist[4]
			amino_acid = linelist[5]
			if gene_name in amino_acid_dict[clone_name]:
				amino_acid_dict[clone_name][gene_name].append(amino_acid)
				mutation_dict[clone_name][gene_name].append(mutation)
			else:
				amino_acid_dict[clone_name][gene_name] = [amino_acid]
				mutation_dict[clone_name][gene_name] = [mutation]

##Second loop to identify any amino acids with two mutations in the same pop
filename_out = 'mutation_data/mutation_lists_with_aa_positions_reannotated.txt'
file_out = open(filename_out,'w')
for line in file_lines:
	correction_tag = 0
	linelist = line.strip().split('\t')
	if len(linelist) < 1.5:
		#Go to the next clone
		pop = linelist[0]
		
	else:
		mutation = ('_').join(str(i) for i in linelist)
		if len(linelist) > 5.5:
			
			gene_name = linelist[4]
			amino_acid = linelist[5]
			same_aa_inds = numpy.where(numpy.array(amino_acid_dict[pop][gene_name]) == amino_acid)[0]
			
			if len(same_aa_inds) > 1.5:
				correction_tag = 1
				
				mut1_list = mutation_dict[pop][gene_name][same_aa_inds[0]].split('_')
				mut2_list = mutation_dict[pop][gene_name][same_aa_inds[1]].split('_')
				
				chr = mut1_list[0]
				loc1 = int(mut1_list[1])
				loc2 = int(mut2_list[1])
				ref1 = mut1_list[2]
				ref2 = mut2_list[2]
				substitution1 = mut1_list[3].strip()
				substitution2 = mut2_list[3].strip()
				
				##Re-translate the amino acid
				
				flanking_region = sequence[chr][0][loc1-11:loc1+10]
				
				if strand[gene_name] == '+':
					loc1_in_codon = (loc1 - gene_start_stop[gene_name][0]- float(phase[gene_name]))%3
					
					loc1_in_codon = int(loc1_in_codon)
					codon_number = numpy.ceil((loc1 - gene_start_stop[gene_name][0]- float(phase[gene_name]) + 1)/3.)
					codon = flanking_region[10 - loc1_in_codon:13 - loc1_in_codon]
					print codon
					print substitution1
					print substitution2
					relative_pos = loc2 - loc1
					print relative_pos
					amino_acid = genetic_code[codon]
					
					if relative_pos == 1:
					
						new_codon = ('').join((codon[0:loc1_in_codon],substitution1,substitution2,codon[loc1_in_codon+1+relative_pos:]))
					else:
						new_codon = ('').join((substitution1,codon[1],substitution2))

					
					new_amino_acid = genetic_code[new_codon]
					
					print '+ strand'
					print gene_name
					print codon_number
					print amino_acid, new_amino_acid
					
				elif strand[gene_name] == '-':
			
					loc1_in_codon = (gene_start_stop[gene_name][1]- loc1 - float(phase[gene_name]))%3
					loc1_in_codon = int(loc1_in_codon)
					codon_number = numpy.ceil((gene_start_stop[gene_name][1] - loc1 - float(phase[gene_name])+ 1)/3.)
					codon = flanking_region[10-(2-loc1_in_codon):13-(2-loc1_in_codon)][::-1]
					ccodon = string.translate(codon,complem)
					relative_pos = loc2 - loc1
					amino_acid = genetic_code[ccodon]
					print codon
					print ccodon
					print substitution1
					print substitution2
					relative_pos = loc2 - loc1
					print relative_pos
					
					if relative_pos == 1:
						new_codon = ('').join((codon[0:loc1_in_codon - relative_pos], substitution2, substitution1, codon[loc1_in_codon+1:]))
					else:
						new_codon = ('').join((substitution2,codon[1],substitution1))

					new_ccodon = string.translate(new_codon,complem)
					new_amino_acid = genetic_code[new_ccodon]
					print '- strand'
					print gene_name
					print codon_number
					print amino_acid, new_amino_acid	
			
	if correction_tag == 0:
		file_out.write(line)
	else:
		if len(linelist) > 7.5:
			new_linelist = linelist[0:7]
			print linelist
			freq = linelist[-1].split(',')[1]
			new_entry = ('').join((amino_acid,'-->',new_amino_acid, ',' ,freq))
			new_linelist.append(new_entry)
			file_out.write(('\t').join(new_linelist) + '\n')
		else:
			new_linelist = list(linelist)
			freq = linelist[6].split(',')[1]
			new_linelist[6] = 'Non'
			new_entry = ('').join((amino_acid,'-->',new_amino_acid, ',' ,freq))
			new_linelist.append(new_entry)
			file_out.write(('\t').join(new_linelist) + '\n')
file_out.close()