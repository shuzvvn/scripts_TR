#!/usr/bin/Rscript

# Usage: 
# run_DECIPHER.R \
# --infile=/mnt/c/Users/stc/project/TR02/source_data/aa.fasta/1A12A.fasta \
# --outfile=/mnt/c/Users/stc/project/TR02/TR02.09/1A12A.tsv

# load libraries

library(optparse)

# specify options in a list
option_list = list(
	make_option(c("-i", "--infile"), type="character", default="NA", help="input seq (Required)", metavar="filename"),
	make_option(c("-o", "--outfile"), type="character", default="NA", help="output tsv (Required)", metavar="filename"),
	make_option(c("-n", "--nt"), action="store_true", default=FALSE, help="input is DNA seq [default=%default]")
); 

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(usage = "%prog [options] [--infile=FILENAME] [--outfile=FILENAME]", option_list=option_list)) 


######################## main ########################
library(DECIPHER)

# read input seq file
if ( opt$nt ) {
	seq <- readDNAStringSet(opt$infile)
} else {
	seq <- readAAStringSet(opt$infile)
}

# print input filename
message <- cat("Input:", opt$infile)
write(message, stderr())

message <- cat("length:", length(seq))
write(message, stderr())

# run DetectRepeats
result <- DetectRepeats(seq, processors=NULL, minScore = 0)

# write output
write.table(result, file=opt$outfile, quote=FALSE, sep='\t')
