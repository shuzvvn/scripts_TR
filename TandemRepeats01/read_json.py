import json
import os
import sys
from Bio import SeqIO

path = '/home/stc/Documents/project/annotation/'
in_fasta = '/home/stc/Documents/project/RepeatsDB_36.fasta'

record_dict = SeqIO.to_dict(SeqIO.parse(in_fasta, "fasta"))

print("ID\tseq_len\tstart\tend")
for file in os.listdir(path):
	if file.endswith(".json"):
		current = os.path.join(path, file)
		ID = file.split(".")[0]
		# Opening JSON file
		f = open(current,)
		# returns JSON object as
		# a dictionary
		data = json.load(f)
		# Iterating through the json list
		dic_repeat = data[0]['repeatsdb_consensus_one'][0]
		seq_len = len(record_dict[ID])
		print("%s\t%i\t%i\t%i" %(ID, seq_len, dic_repeat['start'], dic_repeat['end']))
		# Closing file
		f.close()