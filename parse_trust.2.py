import os
import sys
import re


path_predict = '/Users/stc/project/TR02/TR02.06/TRUST.aa/'
out_file = "/Users/stc/project/TR02/TR02.06/TRUST/aa.tsv"

out_f = open(out_file, "w")
out_f.write("seqID\tstart\tend\tscore1\tscore2\n")



# get predict
for file in os.listdir(path_predict):
	if file.endswith(".out"):
		current = os.path.join(path_predict, file)
		seqID = file.split(".")[0]
		# Opening txt file
		with open(current) as f:
			contents = f.read()
			list_info = re.findall(r'([0-9]+)..([0-9]+) scored ([\S]+) scored \(no gaps\) ([\S]+) included', contents)
			for repeat in list_info:
				start = int(repeat[0])
				end = int(repeat[1])
				score1 = float(repeat[2])
				score2 = float(repeat[3])
				out_f.write('%s\t%i\t%i\t%f\t%f\n' % (seqID, start, end, score1, score2))
				print('%s\t%i\t%i\t%f\t%f' % (seqID, start, end, score1, score2))
out_f.close()