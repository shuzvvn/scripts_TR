import os
import sys
import re

path_ID = '/mnt/c/Users/vivia/EWLab/ID.1511.txt'

with open(path_ID) as f:
    list_ID = f.read().splitlines()

out_f = open("/mnt/c/Users/vivia/EWLab/trust/trust.tsv", "w")
out_f.write("ID\tIndex\tBegin\tEnd\tScore1\tScore2\n")

path_predict = '/mnt/c/Users/vivia/EWLab/trust/BLOSUM62/'


# get predict
for file in os.listdir(path_predict):
	if file.endswith(".out"):
		current = os.path.join(path_predict, file)
		ID = file.split(".")[0]
		Index = list_ID.index(ID) + 1
		# Opening txt file
		with open(current) as f:
			contents = f.read()
			list_info = re.findall(r'([0-9]+)..([0-9]+) scored ([\S]+) scored \(no gaps\) ([\S]+) included', contents)
			for repeat in list_info:
				start = int(repeat[0])
				end = int(repeat[1])
				score1 = float(repeat[2])
				score2 = float(repeat[3])
				out_f.write('%s\t%i\t%i\t%i\t%f\t%f\n' % (ID, Index, start, end, score1, score2))
				print('%s\t%i\t%i\t%i\t%f\t%f' % (ID, Index, start, end, score1, score2))
out_f.close()