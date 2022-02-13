import os
import sys
import re
import pandas as pd
import statistics


# evaluate True Positive
# number of correctly predicted residues / total number of residues in true repeat regions
def evaluate_TP(seq_len, answer_regions, predict_regions):
	# store infor as a list of Y and N
	answer = ['N']*seq_len
	predict = ['N']*seq_len

	# Y: this pos/aa is in the repeat region
	# N: not in repeat region
	for answer_region in answer_regions:
		answer_range = range(answer_region[0]-1, answer_region[1])
		for pos in answer_range:
			answer[pos] = 'Y'

	for predict_region in predict_regions:
		predict_range = range(int(predict_region[0])-1, int(predict_region[1]))
		for pos in predict_range:
			predict[pos] = 'Y'

	# get TP %
	TP_numer = 0
	TP_denom = 0
	for pos in range(seq_len):
		if answer[pos]=='Y':
			TP_denom += 1
			if predict[pos]=='Y':
				TP_numer += 1

	TP_rate = TP_numer/TP_denom

	return TP_rate


# evaluate False Positive
# number of wrongly predicted residues / total number of residues in decoy set
def evaluate_FP(seq_len, predict_regions):
	FP_numer = 0
	FP_denom = seq_len
	for predict_region in predict_regions:
		FP_numer += (int(predict_region[1])-int(predict_region[0])+1)

	FP_rate = FP_numer/FP_denom

	return FP_rate


# get predict
def get_dict_predict(path_predict, index_seqID):
	dict_predict = {}
	df = pd.read_csv(path_predict, sep='\t')
	with open(index_seqID) as f_index_seqID:
		lines = f_index_seqID.readlines()
		for i in range(len(lines)):
			ID_in_decipher = i+1
			seqID = lines[i].rstrip('\n').split('\t')[0]
			df_h = df[(df['Index'] == ID_in_decipher)]
			starts = df_h['Begin'].tolist()
			ends = df_h['End'].tolist()
			list_pos = []
			for region_h in range(len(starts)):
				list_pos.append(tuple([starts[region_h],ends[region_h]]))
			dict_predict[seqID] = list_pos
	return dict_predict


# get answer
def get_dict_answer(path_answer, index_seqID):
	dict_answer = {}
	df = pd.read_csv(path_answer, sep='\t')
	with open(index_seqID) as f_index_seqID:
		lines = f_index_seqID.readlines()
		for i in range(len(lines)):
			seqID = lines[i].rstrip('\n').split('\t')[0]
			df_h = df[(df['ID1'] == seqID)]
			starts = df_h['begin_label'].tolist()
			ends = df_h['end_label'].tolist()
			list_pos = []
			for region_h in range(len(starts)):
				list_pos.append(tuple([starts[region_h],ends[region_h]]))
			dict_answer[seqID] = list_pos
	return dict_answer


dir_aa = '/mnt/c/Users/stc/project/DECIPHER.aa/'
dir_decoy = '/mnt/c/Users/stc/project/DECIPHER.decoy/'

path_answer = '/mnt/c/Users/stc/project/convert_seq_id.tsv'
index_seqID = '/mnt/c/Users/stc/project/seq_len.tsv'

out_file = '/mnt/c/Users/stc/project/decipher_evaluate.1.tsv'
with open(out_file, 'w') as out_file_h:
	out_file_h.write("cutoff\tTPR\tFPR\n")
	print("cutoff\tTPR\tFPR")
	for filename in os.listdir(dir_aa):
		path_aa = dir_aa + filename
		path_decoy = dir_decoy + filename


		aa_dict_predict = get_dict_predict(path_aa, index_seqID)
		decoy_dict_predict = get_dict_predict(path_decoy, index_seqID)

		dict_answer = get_dict_answer(path_answer, index_seqID)


		# evaluate, get true and false positive rates
		TPR_l = []
		FPR_l = []
		with open(index_seqID) as f_index_seqID:
			lines = f_index_seqID.readlines()
			for i in range(len(lines)):
				words = lines[i].rstrip('\n').split('\t')
				seqID = words[0]
				seq_len = int(words[1])

				answer_regions = dict_answer[seqID]
				aa_predict_regions = aa_dict_predict[seqID]
				decoy_predict_regions = decoy_dict_predict[seqID]
				
				TPR = evaluate_TP(seq_len, answer_regions, aa_predict_regions)
				FPR = evaluate_FP(seq_len, decoy_predict_regions)

				#print('%s\t%i\t%.4f\t%.4f' % (seqID, seq_len, TPR, FPR))
				#out_file_h.write('%s\t%i\t%i\t%i\t%i\t%i\t%.4f\t%.4f\n' % (seqID, seq_len, values_d["TP"], values_d["FP"], values_d["TN"], values_d["FN"], TPR, FPR ))
				TPR_l.append(TPR)
				FPR_l.append(FPR)
		TPR_mean = statistics.mean(TPR_l)
		FPR_mean = statistics.mean(FPR_l)
		print('%s\t%.4f\t%.4f' %(filename, TPR_mean, FPR_mean))
		out_file_h.write('%s\t%.4f\t%.4f\n' %(filename, TPR_mean, FPR_mean))