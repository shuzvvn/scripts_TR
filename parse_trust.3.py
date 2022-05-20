#!/usr/bin/env python3

# parse_trust.3.py

# Shu-Ting Cho <vivianlily6@hotmail.com>
# Read TRUST's result and output tsv

# v1 2022/04/21
# v2 2022/03/02
# v3 2022/05/17 Use argparse


# Usage:
# python3 parse_trust.3.py --in_file=~/EW_project/TR02.24/run1/aa_TRUST/K24263.out --out_file=~/EW_project/TR02.24/run1/aa_tsv/K24263.tsv

## import modules
import os
import sys
import re
import argparse

# parsing arguments
def parse_args():
	parser = argparse.ArgumentParser(description="Read TRUST's result and output tsv.")
	parser.add_argument('-i', '--in_file', help="TRUST's result")
	parser.add_argument('-o', '--out_file', help="output .tsv file")
	return parser.parse_args()

args = parse_args()

out_f = open(str(args.out_file), "w")
out_f.write("seqID\tstart\tend\tscore1\tscore2\n")

# Opening result file
with open(args.in_file) as in_f:
	lines = in_f.readlines()
	for line in lines:
		if line[0:2] == '>K':
			line = line.rstrip()
			seqID = line[1:]
			
		elif line[:6] == '# For ' and line[-9:-1] == 'included':
			line = line.rstrip()
			list_info = re.findall(r'([0-9]+)..([0-9]+) scored ([\S]+) scored \(no gaps\) ([\S]+) included', line)
			print(list_info[0])
			start = int(list_info[0][0])
			end = int(list_info[0][1])
			score1 = float(list_info[0][2])
			score2 = float(list_info[0][3])
			out_f.write('%s\t%i\t%i\t%f\t%f\n' % (seqID, start, end, score1, score2))
out_f.close()