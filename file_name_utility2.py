import sys

def sample_list():
	file_in = 'sample_to_filenames_cluster.txt'
	file = open(file_in, 'r')

	unique_file_roots = []
	breseq_input_sample_paths = {}

	firstline = True
	for line in file:
		if firstline:
			firstline = False
		else:
			linelist = line.strip().split('\t')
			file_paths = linelist[-1].split(',')
			sample = linelist[0]
			breseq_input_sample_paths[sample] = []
		
			for file in file_paths:
			
				file_name = file.split('/')[-1]
				file_name_pieces = file_name.split('.')
				file_name_root = file_name.split('.')[0]
			
				breseq_input_sample_paths[sample].append( 'trimmed_fastq_files/' + file_name_root + '.R1.paired' )
				breseq_input_sample_paths[sample].append( 'trimmed_fastq_files/' + file_name_root + '.R2.paired' )
				breseq_input_sample_paths[sample].append( 'trimmed_fastq_files/' + file_name_root + '.R1.unpaired' )
				breseq_input_sample_paths[sample].append( 'trimmed_fastq_files/' + file_name_root + '.R2.unpaired' )
				
				if file_name_root not in unique_file_roots:
					unique_file_roots.append(file_name_root)
		
			breseq_input_sample_paths[sample] = list(set(breseq_input_sample_paths[sample]))
			
	return sorted(breseq_input_sample_paths.keys())