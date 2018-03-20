import numpy
import sys
import matplotlib.pylab as pt

import random
import string
complem = string.maketrans('GATC','CTAG')

#This code takes a .gff file, a .fasta file, and a tab-separated text file with mutations, and outputs a tab-separated text file where the mutations have been annotated.
#Edited 11/12/2015 to handle the case of putative proteins, which lack a 'common' gene name, so that the systematic name will be used.

#read in the genetic code dictionary
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

sequence = []

read_list = []

gene_list = []
gene_start_stop = []
gene_length = []
strand = []
phase = []
gene_list_file = 'w303_ref_with_seq.gff'

file_genes = open(gene_list_file,'r')

file_line_list = list(file_genes)

for i in range(len(file_line_list)):
	if file_line_list[i].startswith('#'):
		continue
	elif file_line_list[i].startswith('>'):
		break
	else:
		#linelist_prev = file_line_list[i-1].split('\t')
		linelist = file_line_list[i].split('\t')
		#print linelist
		if linelist[2] == 'gene':
			#record the gene name, whether it is 
			gene_name = linelist[8].strip().split(';')[1].split('=')[1]
			entry_label = linelist[8].strip().split(';')[1].split('=')[0]
			
		if linelist[2] == 'CDS':
			#record the location and phase of particular CDS's within the gene
			gene_start_stop.append([linelist[0],float(linelist[3]),float(linelist[4])])
			gene_length.append(float(linelist[4]) - float(linelist[3]))
			parent = linelist[8].strip().split(';')[1].split('=')[1].split('_')[0]
			
			if parent == gene_name:
				if entry_label != 'gene':
					gene_list.append(parent)
				else:
					gene_list.append(gene_name)
			else:
				gene_list.append(linelist[8].strip().split(';')[2].split('=')[1])
			
				
			strand.append(linelist[6])
			phase.append(linelist[7])

num_features = len(gene_start_stop)

genome_location = 0

fac = -1

for j in numpy.arange(i,len(file_line_list)):
	if file_line_list[j].startswith('>'):
		fac += 1
		firstline = 1
		continue
	if firstline == 1:
		sequence.append(file_line_list[j].strip())
		firstline = 0
	else:
		sequence[fac] = ('').join((sequence[fac],file_line_list[j].strip()))
file_genes.close()

fasta_chrom = chrom_labels[0]
fasta_chrom_counter = 0

stuck = 0
counter = 0
firstline = 1

#file_out = open('sequencing_data\\all_muts_binary_annotated.txt','w')
file_out = open('mutation_data/all_muts_alts_covs_annotated_with_aa_positions.txt','w')
file_in = open('mutation_data/all_muts_alts_covs.txt','r')

for line in file_in:
		
	# if firstline:
# 		file_out.write(line)
# 		firstline = 0
		
	if line.startswith('chrMito'):
		break
	else:	
		line_list = line.strip().split(',')	
		read_loc = int(line_list[1])
		chrom = line_list[0]
		
		loc_past = 0
		in_a_gene = 0
		substitution = line_list[3].strip()
		anc_allele = line_list[2].strip()
		
		while loc_past < .5:
			if len(gene_start_stop) > genome_location:
				if numpy.logical_and(gene_start_stop[genome_location][1] <= read_loc,gene_start_stop[genome_location][2] >= read_loc):					
					in_a_gene = 1
					loc_past += 1
				elif numpy.logical_and(read_loc < gene_start_stop[genome_location][1],chrom == gene_start_stop[genome_location][0]):
					loc_past += 1
				elif numpy.logical_and(read_loc > gene_start_stop[genome_location][2] ,chrom == gene_start_stop[genome_location][0]):
					genome_location += 1
				elif read_loc < gene_start_stop[genome_location][1]:
					genome_location += 1
				elif read_loc > gene_start_stop[genome_location][2]:
					loc_past += 1
			else:
				loc_past = 1
				
		if in_a_gene < .5:
			file_out.write((line.strip() + '\t\t' + 'Intergenic' + '\t\n'))
		if in_a_gene > .5:
			
			while fasta_chrom != chrom:
					
				fasta_chrom_counter += 1
				fasta_chrom = chrom_labels[fasta_chrom_counter]
			
			flanking_region = sequence[fasta_chrom_counter][read_loc-11:read_loc+10]
			if len(substitution.split(':')) > 1.5:
				fs_indicator = 3
				sub_length=numpy.abs(int(substitution.split(':')[0]))
			else:
				fs_indicator = 0
				sub_length=0
			
			if (strand[genome_location] == '+' and fs_indicator < 2.5):
				loc_in_codon = (read_loc - gene_start_stop[genome_location][1]- float(phase[genome_location]))%3
				loc_in_codon = int(loc_in_codon)
				codon_number = numpy.ceil((read_loc - gene_start_stop[genome_location][1]- float(phase[genome_location]) + 1)/3.)
				codon = flanking_region[10 - loc_in_codon:13 - loc_in_codon]
				
				amino_acid = genetic_code[codon]
				new_codon = ('').join((flanking_region[10-loc_in_codon:10],substitution,flanking_region[11:13-loc_in_codon]))
				
				#new_codon = new_sequence[:3]
				new_amino_acid = genetic_code[new_codon]
			elif strand[genome_location] == '+':
				codon_number = numpy.ceil((read_loc - gene_start_stop[genome_location][1]- float(phase[genome_location]) + 1)/3.)
			elif (strand[genome_location] == '-' and fs_indicator < 2.5):
				
				
				loc_in_codon = (gene_start_stop[genome_location][2]- read_loc - float(phase[genome_location]))%3
				loc_in_codon = int(loc_in_codon)
				codon_number = numpy.ceil((gene_start_stop[genome_location][2] - read_loc - float(phase[genome_location])+ 1)/3.)
				codon = flanking_region[10-(2-loc_in_codon):13-(2-loc_in_codon)][::-1]
				ccodon = string.translate(codon,complem)
				
				amino_acid = genetic_code[ccodon]
				new_codon = ('').join((flanking_region[10-(2-loc_in_codon):10],substitution,flanking_region[11:13-(2-loc_in_codon)]))[::-1]
				#new_codon = new_sequence[:3]
				new_ccodon = string.translate(new_codon,complem)
				new_amino_acid = genetic_code[new_ccodon]
			else:
				codon_number = numpy.ceil((gene_start_stop[genome_location][2] - read_loc - float(phase[genome_location]) + 1)/3.)
				
			if (fs_indicator > 2.5 and sub_length < 2.5):
				file_out.write((line.strip() + '\t' + gene_list[genome_location] + '\t' + str(codon_number) + '\t' + 'Non' + '\t' + 'FS' + '\n'))
			elif (fs_indicator > 2.5 and sub_length == 3):
				file_out.write((line.strip() + '\t' + gene_list[genome_location] + '\t' + str(codon_number) + '\t' +'Non' + '\t' + 'inframe_indel' + '\n'))
			elif new_amino_acid == amino_acid:
				file_out.write((line.strip() + '\t' + gene_list[genome_location] +'\t' + str(codon_number) + '\t' +'Syn' + '\t\n'))
			else:
				file_out.write((line.strip() + '\t' + gene_list[genome_location] + '\t' + str(codon_number) + '\t' + 'Non' + '\t'+ ('').join((amino_acid,'-->',new_amino_acid)) + '\n'))
			
file_in.close()
file_out.close()

		