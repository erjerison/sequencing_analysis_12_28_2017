import sys
import os
import numpy

def get_chrom_dict():
	ref_idx = '~/w303_reference_genome/w303_ref.fasta.fai'

	file=open(ref_idx,'r')
	chrom_dict={}
	for line in file:
		linestrs = line.split('\t')
		length = numpy.float(linestrs[1])
		
		chrom_dict[linestrs[0]] = length

	return chrom_dict

if __name__=='__main__':
	if sys.argv[1]=='keys':
		chrom_dict = get_chrom_dict()
		for key in chrom_dict:
			print key
	else:
		chrom_dict = get_chrom_dict()
		length = int(chrom_dict[sys.argv[1]])
		for j in numpy.arange(0,length+50000,50000):
			print j
