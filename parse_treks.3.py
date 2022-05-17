#!/usr/bin/env python3

# parse_treks.3.py

# Shu-Ting Cho <vivianlily6@hotmail.com>
# Read T-REKS's result and output tsv

# v1 2022/04/21
# v2 2022/03/02
# v3 2022/05/12 Use argparse


# Usage:
# python3 parse_treks.3.py --in_file=~/project/TR02/TR02.23/aa_TREKS/K24263.out --out_file=~/project/TR02/TR02.23/aa_tsv/K24263.tsv

## import modules
import os
import sys
import re
import argparse

# parsing arguments
def parse_args():
	parser = argparse.ArgumentParser(description="Read T-REKS's result and output tsv.")
	parser.add_argument('-i', '--in_file', help="T-REKS's result")
	parser.add_argument('-o', '--out_file', help="output .tsv file")
	return parser.parse_args()


args = parse_args()

out_f = open(str(args.out_file), "w")
out_f.write("seqID\tstart\tend\tscore\n")


# Opening result file
with open(str(args.in_file)) as in_f:
	lines = in_f.readlines()
	for line in lines:
		if line[0] == '>':
			line = line.rstrip()
			seqID = line[1:]
			
		elif line[:6] == 'Length':
			line = line.rstrip()
			list_info = re.findall(r'from[\s]+([0-9]+)[\s]+to[\s]+([0-9]+)[\s]+-[\s]+Psim:([\S]+)', line)
			start = int(list_info[0][0])
			end = int(list_info[0][1])
			score = float(list_info[0][2])
			out_f.write('%s\t%i\t%i\t%f\n' % (seqID, start, end, score))
			#print('%s\t%i\t%i\t%f' % (seqID, start, end, score))
out_f.close()