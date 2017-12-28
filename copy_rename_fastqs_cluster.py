from shutil import copyfile

file_in = open('sample_to_filenames_cluster.txt','r')

firstline = 1

line_counter = 0

for line in file_in:
	
	if firstline:
		firstline = 0
		
	else:
		line_list = line.strip().split('\t')
		old_filenames = line_list[-2].split(',')
		new_filenames = line_list[-1].split(',')
		
		for i in range(len(old_filenames)):
		
			old_filename = old_filenames[i]
			new_filename = new_filenames[i]
			
			copyfile(old_filename, new_filename)
		
		line_counter += 1
		print 'Done with sample', line_counter

file_in.close()