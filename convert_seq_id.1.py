import os
import sys
import re
import pandas as pd
import numpy as np
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

path = '/mnt/c/Users/stc/project/mmcif/'
in_tsv = '/mnt/c/Users/stc/project/cesymm.tsv'

print('ID1\tID2\tseq_len\tbegin_auth\tend_auth\tbegin_label\tend_label')

# get cesymm result
df = pd.read_csv(in_tsv, sep='\t')

for row in range(df.shape[0]):
	seqID1 = df.iloc[row, 0]
	seqID = df.iloc[row, 1]
	seq_len = df.iloc[row, 2]
	result_h = df.iloc[row, 3]
	result_h = result_h.replace(seqID+'_', '') # remove seqID

	# read mmcif
	pdbID, chainID = seqID.split('.')
	file_cif = pdbID + '.cif'
	mmcif_dict = MMCIF2Dict(path+file_cif)

	label_seq_ids = mmcif_dict['_atom_site.label_seq_id']
	auth_seq_ids = mmcif_dict['_atom_site.auth_seq_id']
	auth_asym_ids = mmcif_dict['_atom_site.auth_asym_id'] # chain ID list
	data = np.column_stack([auth_asym_ids, auth_seq_ids, label_seq_ids])
	new_array = [tuple(row) for row in data]
	uniques = np.unique(new_array, axis=0)

	lookup = {}
	for chain_h in set(auth_asym_ids):
		lookup[chain_h] = {}
		array_h = uniques[uniques[:, 0] == chain_h]

		for i in range(array_h.shape[0]):
			lookup[chain_h][array_h[i,1]] = array_h[i,2]

	copies = result_h.split(';')
	for copy in copies:
		copy = re.sub(r'([0-9])-', r'\1@', copy)
		begin_auth, end_auth = copy.split('@')
		begin_label, end_label = lookup[chainID][begin_auth], lookup[chainID][end_auth]
		print('%s\t%s\t%s\t%s\t%s\t%s\t%s' % (seqID1, seqID, seq_len, begin_auth, end_auth, begin_label, end_label))
