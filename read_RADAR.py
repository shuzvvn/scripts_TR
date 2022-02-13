import os
import sys
import re


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




path_predict = 'RADAR/'
path_answer = 'RepeatsDB_36.answer.txt'

# get predict
dict_predict = {}
for file in os.listdir(path_predict):
	if file.endswith(".txt"):
		current = os.path.join(path_predict, file)
		ID = file.split(".")[0]
		# Opening txt file
		with open(current) as f:
			contents = f.read()
			# pattern for start and end: ^ +([0-9]+)- +([0-9]+)
			list_pos = re.findall(r'\s+([0-9]+)-\s+([0-9]+)', contents)
			dict_predict[ID] = list_pos

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
print('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % ('ID', 'seq_len', 'TP', 'FP', 'TN', 'FN','TPR','FPR'))
for ID in dict_answer:
	seq_len = dict_answer[ID][0]
	answer_region = [dict_answer[ID][1], dict_answer[ID][2]]
	predict_regions = dict_predict[ID]
	
	values_d = evaluate(seq_len, answer_region, predict_regions)

	# true positive rate, sensitivity (TPR = TP/(TP+FN))
	TPR = values_d["TP"] / (values_d["TP"]+values_d["FN"])
	# False Positive Rate (FPR), specificity ( = FP/(TN+FP))
	FPR = values_d["FP"] / (values_d["TN"]+values_d["FP"])

	print('%s\t%i\t%i\t%i\t%i\t%i\t%.4f\t%.4f' % (ID, seq_len, values_d["TP"], values_d["FP"], values_d["TN"], values_d["FN"], TPR, FPR ))

