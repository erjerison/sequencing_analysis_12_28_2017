import numpy

N7indices = ['TAAGGCGA','CGTACTAG','AGGCAGAA','TCCTGAGC','GGACTCCT','TAGGCATG','CTCTCTAC','CAGAGAGG','GCTACGCT','CGAGGCTG','AAGAGGCA','GTAGAGGA']
N5indices = ['TAGATCGC', 'CTCTCTAT', 'TATCCTCT', 'AGAGTAGA', 'GTAAGGAG', 'ACTGCATA', 'AAGGAGTA', 'CTAAGCCT']

rows = ['A','B','C','D','E','F','G','H']

N5_dict = dict(zip(rows, N5indices))
N7_dict = dict(zip(numpy.arange(12), N7indices))

file1 = open('plate_layouts_for_sequencing_updated.txt','r')

file_out = open('sample_to_index_key.txt','w')
file_out.write( ('\t').join( ('Clone','index','library pool') ) + '\n' )

for line in file1:
	if line.startswith('#'):
		continue
	elif line.startswith('>'):
		plate = line.strip()[1:]
		continue
	else:
		row_list = line.strip().split('\t')
		row = row_list[0]
		
		N5ind = N5_dict[row]
		
		for j in numpy.arange(len(row_list)-1):
			
			sample = row_list[j+1]
			
			N7ind = N7_dict[j]
			
			index_tag = N7ind + '-' + N5ind
		
			if plate == 'Plate1':
				library_pool = 1
			
			elif (plate == 'Plate2' and j == 9): ##Col 10 from plate 2 was put with plate 3, which only included cols 1-9, to balance the number of samples per lane
				library_pool = 3
					
			elif plate == 'Plate2':
				library_pool = 2
			
			else:
				library_pool = 3
			
			if sample != 'b': #These are intentional blanks
				file_out.write( ('\t').join( (sample, index_tag, str(library_pool)) ) + '\n')

file1.close()
file_out.close()
			