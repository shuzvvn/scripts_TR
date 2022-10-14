#!/usr/bin/env Rscript

# Usage: 
# DECIPHER_DetectRepeats.3.r \
# --infile=~/EW_project/Kgroup/nt/K24316.fas.gz \
# --DetectRepeat=~/EW_project/TR02.16/DetectRepeats_nt/K24316.DetectRepeats.rds \
# --type="tandem"
# --nt

# load libraries

library(optparse)

# specify options in a list
option_list = list(
	make_option(c("-i", "--infile"), type="character", default="NA", help="input seq (Required)", metavar="filename"),
	make_option(c("-n", "--nt"), action="store_true", default=FALSE, help="input is DNA seq [default=%default]"),
	make_option(c("-t", "--type"), type="character", default="tandem", help="Character string indicating the type of repeats to detect. This should be (an abbreviation of) one of 'tandem', 'interspersed', or 'both'. [default=%default]", metavar="string")
); 

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(usage = "%prog [options] [--infile=FILENAME] [--outfile=FILENAME]", option_list=option_list)) 


######################## main ########################
suppressMessages(library(DECIPHER))

source('DetectRepeats.R')
source('ScoreAlignment.R')
source('utils.R')

# read input seq file
if ( opt$nt ) {
	# print input filename
	message <- cat("readDNAStringSet", opt$infile)
	write(message, stderr())
	seq <- readDNAStringSet(opt$infile)
} else {
	# print input filename
	message <- cat("readAAStringSet", opt$infile)
	write(message, stderr())
	seq <- readAAStringSet(opt$infile)
}

# run DetectRepeats
result_DetectRepeats <- DetectRepeats(seq, processors=NULL, verbose=TRUE, type = opt$type)
write('\nDetectRepeats DONE!', stderr())

# write output
saveRDS(result_DetectRepeats, file=opt$outfile, compress = TRUE)