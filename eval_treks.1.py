import os
import pandas as pd
import statistics

# ref info
path_answer = '/Users/stc/project/TR02/source_data/convert_seq_id.tsv'
seqID_len_file = '/Users/stc/project/TR02/source_data/seq_len.tsv'

# detector results
dir_aa = '/Users/stc/project/TR02/TR02.05/treks.eval.aa/'
dir_decoy = '/Users/stc/project/TR02/TR02.05/treks.eval.decoy/'

# output
out_file = '/Users/stc/project/TR02/TR02.05/treks.evaluate.1.tsv'


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


# get list of positions for each seqID as dict
def get_pos_dict(file_path, seqID_len_df, colname_seqID='seqID', colname_start='start', colname_end='end'):
	dict_predict = {}
	df = pd.read_csv(file_path, sep='\t')
	for seqID in seqID_len_df.index.values:
		df_h = df[(df[colname_seqID] == seqID)]
		starts = df_h[colname_start].tolist()
		ends = df_h[colname_end].tolist()
		list_pos = []
		for region_i in range(len(starts)):
			list_pos.append(tuple([starts[region_i],ends[region_i]]))
		dict_predict[seqID] = list_pos
	return dict_predict


with open(out_file, 'w') as out_file_h:
	out_file_h.write("cutoff\tTPR\tFPR\n")
	print("cutoff\tTPR\tFPR")
	for filename in os.listdir(dir_aa):
		path_aa = dir_aa + filename
		path_decoy = dir_decoy + filename

		# get seqID and len info
		seqID_len_df = pd.read_csv(seqID_len_file, sep='\t', names=['len'], header=None, index_col=0)

		# get list of positions for each seqID as dict
		aa_dict_predict = get_pos_dict(path_aa, seqID_len_df, colname_seqID='seqID', colname_start='start', colname_end='end')
		decoy_dict_predict = get_pos_dict(path_decoy, seqID_len_df, colname_seqID='seqID', colname_start='start', colname_end='end')

		dict_answer = get_pos_dict(path_answer, seqID_len_df, colname_seqID='ID1', colname_start='begin_label', colname_end='end_label')

		# evaluate, get true and false positive rates
		TPR_l = []
		FPR_l = []
		for seqID in seqID_len_df.index.values: # for each seq
			seq_len = seqID_len_df.loc[seqID,'len']

			# get the predicted and true repeat region positions of this seq
			answer_regions = dict_answer[seqID]
			aa_predict_regions = aa_dict_predict[seqID]
			decoy_predict_regions = decoy_dict_predict[seqID]
			
			TPR = evaluate_TP(seq_len, answer_regions, aa_predict_regions)
			FPR = evaluate_FP(seq_len, decoy_predict_regions)

			TPR_l.append(TPR)
			FPR_l.append(FPR)
		TPR_mean = statistics.mean(TPR_l)
		FPR_mean = statistics.mean(FPR_l)
		print('%s\t%.4f\t%.4f' %(filename, TPR_mean, FPR_mean))
		out_file_h.write('%s\t%.4f\t%.4f\n' %(filename, TPR_mean, FPR_mean))