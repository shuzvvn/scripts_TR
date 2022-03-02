import os
import sys
import re

path_predict = '/mnt/c/Users/stc/project/TR02.05/T-REKS/DECOY.txt'
out_file = '/mnt/c/Users/stc/project/TR02.05/T-REKS.decoy.tsv'

out_f = open(out_file, "w")
out_f.write("seqID\tstart\tend\tscore\n")


# Opening result file
with open(path_predict) as in_f:
	lines = in_f.readlines()
	for line in lines:
		if line[0] == '>':
			line = line.rstrip()
			seqID = line[1:]
			
		elif line[:6] == 'Length':
			line = line.rstrip()
			list_info = re.findall(r'from[\s]+([0-9]+)[\s]+to[\s]+([0-9]+)[\s]+-[\s]+Psim:([\S]+)', line)
			start = int(list_info[0][0])
			end = int(list_info[0][1])
			score = float(list_info[0][2])
			out_f.write('%s\t%i\t%i\t%f\n' % (seqID, start, end, score))
			print('%s\t%i\t%i\t%f' % (seqID, start, end, score))
out_f.close()