import os
import sys
import re

path_ID = '/mnt/c/Users/vivia/EWLab/ID.1511.txt'

with open(path_ID) as f:
    list_ID = f.read().splitlines()

out_f = open("/mnt/c/Users/vivia/EWLab/treks/treks.tsv", "w")
out_f.write("Index\tID\tBegin\tEnd\tScore\n")

path_predict = '/mnt/c/Users/vivia/EWLab/treks/treks.out.txt'

# Opening result file
with open(path_predict) as in_f:
	contents = in_f.readlines()
	for line in contents:
		if line[0] == '>':
			line = line.rstrip()
			ID = line[1:]
			Index = list_ID.index(ID) + 1
		elif line[:6] == 'Length':
			line = line.rstrip()
			list_info = re.findall(r'from[\s]+([0-9]+)[\s]+to[\s]+([0-9]+)[\s]+-[\s]+Psim:([\S]+)', line)
			start = int(list_info[0][0])
			end = int(list_info[0][1])
			score = float(list_info[0][2])
			out_f.write('%s\t%i\t%i\t%i\t%f\n' % (ID,Index, start, end, score))
			print('%s\t%i\t%i\t%i\t%f' % (ID,Index, start, end, score))
out_f.close()