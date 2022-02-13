import os
import sys
import re
import pandas as pd
import statistics


## evaluate
# make a list with seq_length of elements
def evaluate(seq_len, answer_region, predict_regions):
	# store infor as a list of Y and N
	answer = ['N']*seq_len
	predict = ['N']*seq_len

	# Y: this pos/aa is in the repeat region
	# N: not in repeat region
	answer_range = range(answer_region[0]-1, answer_region[1])
	for pos in answer_range:
		answer[pos] = 'Y'

	for predict_region in predict_regions:
		predict_range = range(int(predict_region[0])-1, int(predict_region[1]))
		for pos in predict_range:
			predict[pos] = 'Y'

	# compare
	values_d = {'TP': 0, 'FP': 0, 'TN': 0, 'FN': 0}
	for pos in range(seq_len):
		if predict[pos]=='Y' and answer[pos]=='Y':
			values_d['TP'] += 1
		elif predict[pos]=='N' and answer[pos]=='N':
			values_d['TN'] += 1
		elif predict[pos]=='Y' and answer[pos]=='N':
			values_d['FP'] += 1
		else:
			values_d['FN'] += 1
	return values_d


path_predict = 'decipher/DetectRepeats_1511_v2.tsv'
path_answer = 'RepeatsDB_1511.answer.tsv'
index_seqID = 'ID.1511.txt'
out_file = 'decipher/decipher_evaluate.tsv'


# get predict
dict_predict = {}
df = pd.read_csv(path_predict, sep='\t')
with open(index_seqID) as f_index_seqID:
	lines = f_index_seqID.readlines()
	for i in range(len(lines)):
		ID_in_decipher = i+1
		seqID = lines[i].rstrip('\n')
		df_h = df[(df['Index'] == ID_in_decipher)]
		starts = df_h['Begin'].tolist()
		ends = df_h['End'].tolist()
		list_pos = []
		for region_h in range(len(starts)):
			list_pos.append(tuple([starts[region_h],ends[region_h]]))
		dict_predict[seqID] = list_pos


# get answer
dict_answer = {}
with open(path_answer) as f:
	line = f.readline()
	while line:
		line = f.readline()
		words = line.rstrip('\n').split("\t")
		if len(words)>1:
			# ID, seq_len, start, end
			dict_answer[words[0]] = [int(words[1]), int(words[2]), int(words[3])]


# evaluate, get counts of true and false positive and negative
TPR_l = []
FPR_l = []
with open(out_file, 'w') as out_file_h:
	out_file_h.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % ('ID', 'seq_len', 'TP', 'FP', 'TN', 'FN','TPR','FPR'))
	for ID in dict_answer:
		seq_len = dict_answer[ID][0]
		answer_region = [dict_answer[ID][1], dict_answer[ID][2]]
		predict_regions = dict_predict[ID]
		
		values_d = evaluate(seq_len, answer_region, predict_regions)

		# true positive rate, sensitivity (TPR = TP/(TP+FN))
		TPR = values_d["TP"] / (values_d["TP"]+values_d["FN"])
		# False Positive Rate (FPR), specificity ( = FP/(TN+FP))
		FPR = values_d["FP"] / (values_d["TN"]+values_d["FP"])

		#print('%s\t%i\t%i\t%i\t%i\t%i\t%.4f\t%.4f' % (ID, seq_len, values_d["TP"], values_d["FP"], values_d["TN"], values_d["FN"], TPR, FPR ))
		out_file_h.write('%s\t%i\t%i\t%i\t%i\t%i\t%.4f\t%.4f\n' % (ID, seq_len, values_d["TP"], values_d["FP"], values_d["TN"], values_d["FN"], TPR, FPR ))
		TPR_l.append(TPR)
		FPR_l.append(FPR)

print('TPR = %.4f, FPR = %.4f' %(statistics.mean(TPR_l), statistics.mean(FPR_l)))