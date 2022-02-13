#!/usr/bin/python3

# read_repeatDB_json.2.py

# Shu-Ting Cho <shutingc@andrew.cmu.edu>
# read the repeat region from json

# v1 2021/10/20


# modules
import argparse
import json
import os
import sys
from Bio import SeqIO

# parsing arguments
def parse_args():
	parser = argparse.ArgumentParser(description='Parse the repeat region from json.')
	parser.add_argument('-i', '--in_dir', help='input directory that contains .json files')
	parser.add_argument('-f', '--fasta', help='input fasta files')
	parser.add_argument('-o', '--out_file', default="./out.tsv", help='output tsv (default: current working directory)')
	return parser.parse_args()

# make dir if not exist
def make_outdir(out_file):
	out_dir = os.path.dirname(out_file)
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)


# main
def main():
	args = parse_args()

	path = args.in_dir
	in_fasta = args.fasta
	out_file = args.out_file
	make_outdir(out_file)

	record_dict = SeqIO.to_dict(SeqIO.parse(in_fasta, "fasta"))

	print("ID\tseq_len\tstart\tend")
	with open(out_file, 'w') as out_file_h:
		out_file_h.write("ID\tseq_len\tstart\tend\n")
		for file in os.listdir(path):
			if file.endswith(".json"):
				current = os.path.join(path, file)
				ID = file.split(".")[0]
				# Opening JSON file
				json_f = open(current,)
				# returns JSON object as a dictionary
				data = json.load(json_f)
				# Iterating through the json list
				dic_repeat = data[0]['repeatsdb_consensus_one'][0]
				seq_len = len(record_dict[ID])
				print("%s\t%i\t%i\t%i" %(ID, seq_len, dic_repeat['start'], dic_repeat['end']))
				out_file_h.write("%s\t%i\t%i\t%i\n" %(ID, seq_len, dic_repeat['start'], dic_repeat['end']))

				pdb_chains = data[0]['pdb_chains']
				for record_h in pdb_chains:
					id_h = record_h["id"]
					start_h = record_h["start"]
					end_h = record_h["end"]






				# Closing file
				json_f.close()

if __name__ == '__main__':
	main()