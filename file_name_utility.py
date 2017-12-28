import sys

file_in = 'scripts/sample_to_filenames_cluster.txt'
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
			
			breseq_input_sample_paths[sample].append( 'merged_trimmed_fastqs/' + file_name_root + '.merged.fastq' )
			breseq_input_sample_paths[sample].append( 'merged_trimmed_fastqs/' + file_name_root + '.unmerged_R1.fastq' )
			breseq_input_sample_paths[sample].append( 'merged_trimmed_fastqs/' + file_name_root + '.unmerged_R2.fastq' )
			breseq_input_sample_paths[sample].append( 'trimmed_fastq_files/' + file_name_root + '.R1.unpaired.fastq' )
			breseq_input_sample_paths[sample].append( 'trimmed_fastq_files/' + file_name_root + '.R2.unpaired.fastq' )
				
			if file_name_root not in unique_file_roots:
				unique_file_roots.append(file_name_root)
		
		breseq_input_sample_paths[sample] = list(set(breseq_input_sample_paths[sample]))

if sys.argv[1] == 'list':
	for file_root in unique_file_roots:
		print file_root

elif sys.argv[1] == 'sample_list':
	for sample in sorted(breseq_input_sample_paths.keys()):
		print sample
			
elif sys.argv[1] == 'trial': ##Default behavior: try one job
	print unique_file_roots[0]

else:
	key = sys.argv[1]
	for file_name in breseq_input_sample_paths[key]:
		print file_name