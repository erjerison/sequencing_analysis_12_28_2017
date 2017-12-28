import os, glob
import numpy

##Get the names of all the fastq files

filename_base = '/Volumes/EJ_10817/Desai_lab_data/pleiotropy_project/pleiotropy_sequencing/raw_fastq_files'
new_filename_base = filename_base + '/fastqs_with_sample_names'

run1_lane1_filenames = os.listdir( filename_base + '/plate1_copy/lane1')
run1_lane2_filenames = os.listdir( filename_base + '/plate1_copy/lane2')

run2_foldernames = os.listdir( filename_base + '/plate2_plate3_copy')

run2_lane1_filenames = []
run2_lane2_filenames = []

for folder in run2_foldernames:
	if '7693' in folder:
		path = filename_base + '/plate2_plate3_copy/' + folder
		for file in glob.glob(path + '/*.bz2'):
			run2_lane1_filenames.append( file )
	elif '7694' in folder:
		path = filename_base + '/plate2_plate3_copy/' + folder
		for file in glob.glob(path + '/*.bz2'):
			run2_lane2_filenames.append( file )

##Match sequencing index, run with sample name

file_in = open('sample_to_index_key.txt','r')
file_out = open('sample_to_filenames.txt','w')

file_out.write( ('\t').join(('clone', 'index', 'library pool', 'old filenames', 'new filenames')) + '\n')

firstline = 1

for line in file_in:

	if firstline:
		firstline = 0
		
	else:
		
		old_filenames = []
		new_filenames = []
		
		linelist = line.strip().split()
		sample = linelist[0]
		index = linelist[1]
		library_pool = linelist[2]
	
		if library_pool == '1':
			for name in run1_lane1_filenames:
				if index in name:
					old_filenames.append(filename_base + '/plate1_copy/lane1/' + name)
					name_end = name.split('_')[-1]
					new_filenames.append(filename_base + '/fastqs_with_sample_names/' + ('_').join( (sample, 'run1', 'lane1',name_end) ))
			for name in run1_lane2_filenames:
				if index in name:
					old_filenames.append(filename_base + '/plate1_copy/lane2/' + name)
					name_end = name.split('_')[-1]
					new_filenames.append(filename_base + '/fastqs_with_sample_names/' + ('_').join( (sample, 'run1', 'lane2',name_end) ))
		elif library_pool == '2':
			for name in run2_lane1_filenames:
				if index in name:
					old_filenames.append( name)
					name_end = name.split('_')[-1]
					new_filenames.append(filename_base + '/fastqs_with_sample_names/' + ('_').join( (sample, 'run2', 'lane1', index, name_end) ))
		else:
			for name in run2_lane2_filenames:
				if index in name:
					old_filenames.append( name)
					name_end = name.split('_')[-1]
					new_filenames.append(filename_base + '/fastqs_with_sample_names/' + ('_').join( (sample, 'run2', 'lane2', index, name_end) ))
				print name
				print old_filenames
		file_out.write( ('\t').join( (sample, index, library_pool, (',').join(old_filenames),  (',').join(new_filenames)) ) + '\n')
		
file_in.close()
file_out.close()