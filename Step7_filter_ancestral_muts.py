import numpy
import sys
import matplotlib.pylab as pt
import python_name_change_utility

samples = python_name_change_utility.print_all_names_with_indexes()
	
seg_index_dict = {}
	
for i in range(len(samples)):
	seg_id = samples[i].split('_')[0]
		
	if seg_id in seg_index_dict:
		seg_index_dict[seg_id].append(i)
	else:
		seg_index_dict[seg_id] = []
		seg_index_dict[seg_id].append(i)

print seg_index_dict

chr_trans_dict = {'chrI':'chr01','chrII':'chr02','chrIII':'chr03','chrIV':'chr04','chrV':'chr05','chrVI':'chr06','chrVII':'chr07','chrVIII':'chr08','chrIX':'chr09','chrX':'chr10','chrXI':'chr11','chrXII':'chr12','chrXIII':'chr13','chrXIV':'chr14','chrXV':'chr15','chrXVI':'chr16'}

minority_allele_list = []
n_samples = len(samples)

for chr in chr_trans_dict:
	file_in = 'full_chromosomes/combined_muts_' + chr + '.txt'
	
	file_out_name = 'mutation_data/' + chr_trans_dict[chr] + '_potential_de_novo_muts.txt'
	file = open(file_in, 'r')
	
	file_out = open(file_out_name, 'w')
	file_out2_name = 'mutation_data/' + chr_trans_dict[chr] + '_binary_potential_de_novo_muts.txt'
	
	file_out2 = open(file_out2_name, 'w')
	
	file_out.write('Chr,' + 'Loc,' + 'ref,' + 'alt,' + (' ').join(samples) + ',' + (' ').join(samples) + '\n')
	file_out2.write('Chr,' + 'Loc,' + 'ref,' + 'alt,' + (' ').join(samples) + ',' + '\n')
	
	for line in file:
		entries = line.strip().split(',')
		
		loc_ref_alt = entries[0:4]
		alts = entries[4]
		covs = entries[5]
		
		alt_list = numpy.array([float(a) for a in alts.strip().split()])
		cov_list = numpy.array([float(c) for c in covs.strip().split()])
		
		de_novo = 0
		alt_fracs = alt_list/cov_list
		
		#print alt_fracs
		if (numpy.sum(alt_fracs[~numpy.isnan(alt_fracs)] > .5) > .5 and numpy.nansum(numpy.round(alt_fracs)) < 224.5):
			#first filter: the mutation cannot be at > 50% in 225 or more lines (90% of all lines). This subtracts all ancestral mutations.
			#second filter: filter out alignment artifacts. These come in two types. The first (more common) type is that you sometimes see a minority allele in all the descendants of particular segregants at a particular locus.
			#I believe this happens because there is a nearby indel that is causing alignment errors at a nearby site. This problem may be fixable by creating a new reference for each kruglyak segregant and aligning
			#to that. For now, I am going to filter out all sites where at least 5 different populations show polymorphism in 10% of the reads. This will also take care of the other type of alignment artifact that we usually see,
			#where some regions are just error-prone presumably due to similarity to other regions and hence mis-alignment.
			minority_allele_list.append(numpy.sum(alt_fracs[~numpy.isnan(alt_fracs)] > .1))
			if numpy.sum(alt_fracs[~numpy.isnan(alt_fracs)] > .1) < 4.5:
		
				#print numpy.sum(alt_fracs > .1)
			
				file_out.write(line)
		
				binary_line = ('\t').join(str(a) for a in numpy.round(alt_fracs))
				preamble = (',').join(str(p) for p in loc_ref_alt)
				file_out2.write(preamble + binary_line + '\n')
			
	file_out.close()
	file_out2.close()
	file.close()
print len(minority_allele_list)
pt.hist(minority_allele_list,bins=numpy.arange(0,253))
pt.title('number of pops with allele at >10%')
pt.xlabel('number of alleles')
pt.ylabel('number of pops with this allele at >10%')
pt.xlim([0,253])
pt.show()
				