import numpy
import sys
import file_name_utility

reference_strs = ['.', ',']
snp_strs = ['A','C','T','G', 'N','*'] #, 'a','c','t','g', 'n']
indel_strs = ['+','-']
refskip_strs = ['<','>']
var_strs = []
var_strs.extend(snp_strs)
var_strs.extend(indel_strs)

ref_forward_int = ord('.')
ref_reverse_int = ord(',')
insertion_int = ord('+')
deletion_int = ord('-')

mpileup_file = sys.stdin

#population = sys.argv[1]

# seg_dict = clones_to_wells.clone_to_plate_dict()
# sample_list = []
# for key in sorted(seg_dict.keys()):
# 	for entry in seg_dict[key][1]:
# 		sample_list.append(entry)
# print ", ".join(entry for entry in sample_list)	

sample_list = file_name_utility.print_all_names_with_indexes()
print ", ".join(entry for entry in sample_list)

if len(sys.argv) > 3:
    first_position = long(sys.argv[2])
    last_position = long(sys.argv[3])
else:
    first_position = 0
    last_position = 15000000

num_lines = 0
num_printed_lines = 0

total_mean_depth = 0

for line in mpileup_file:
    num_lines += 1
    if num_lines % 10000 == 0:
        sys.stderr.write("%d lines processed, %d passed\n" % (num_lines,num_printed_lines))

    items = line.split('\t')
    chromosome = items[0]
    position = long(items[1])
    ref = items[2]
    additional_indel_allele = ""
 
    if position < first_position:
        continue
    
    if position > last_position:
        break

    times = []
    alts = []
    depths = []

    num_samples = (len(items)-3)/3

    alleles = {var: [0 for i in xrange(0,num_samples)] for var in var_strs}
        
    for i in xrange(0, num_samples):
       # t = sample_times[i]
            
        depth = long(items[(i+1)*3])
        allele_str = items[(i+1)*3+1].upper()
        depths.append(depth)
       # times.append(t)
        j=0
        while j < len(allele_str):
            if allele_str[j] in reference_strs: # reference
                alleles[ref][i] += 1
            elif allele_str[j] in snp_strs: # regular SNP
                alleles[allele_str[j]][i] += 1
            elif allele_str[j] == '^': # beginning of read
                j+=1 # next character is map quality
            elif allele_str[j] in indel_strs: # indel
                indel_allele = allele_str[j] 
                # correct for counting ref last time
                alleles[ref][i] -= 1        
                j += 1
                k = j
                while allele_str[k].isdigit():
                    k+=1
                indel_len = long(allele_str[j:k])
                j = k
                additional_indel_allele = ('%d:%s' % (indel_len, allele_str[j:j+indel_len]))
                alleles[indel_allele][i] += 1
                j = j + indel_len - 1
            else:
                pass
            j+=1 
    
    depths = numpy.array(depths)
    total_mean_depth += numpy.mean(depths)
     
    for allele in alleles.keys():  
        if allele != ref and allele != '*':
            alts = numpy.array(alleles[allele])

            # alt_map = {}
            # depth_map = {}
            # for t,alt,depth in zip(times,alts,depths):
                # if t not in alt_map:
                    # alt_map[t] = 0
                    # depth_map[t] = 0

                # alt_map[t] += alt
                # depth_map[t] += depth
            
            # merged_times = [t for t in sorted(alt_map.keys())]
            # merged_alts = numpy.array([alt_map[t] for t in merged_times])
            # merged_depths = numpy.array([depth_map[t] for t in merged_times])
            
            if numpy.any((alts >= 4)*(depths>=2)*(alts >= 0.05*depths)):
                if allele in indel_strs:
                    allele = ("%s%s" % (allele,additional_indel_allele))
                print ", ".join([chromosome, str(position), ref, allele, " ".join(str(a) for a in alts), " ".join(str(d) for d in depths)])
	
total_mean_depth = total_mean_depth/float(num_lines)
depth_file = open('avg_depth.txt','a')
depth_file.write(chromosome + ' region ending at ' + str(position) + ',overall mean depth, ' + str(total_mean_depth))
depth_file.close()
	
