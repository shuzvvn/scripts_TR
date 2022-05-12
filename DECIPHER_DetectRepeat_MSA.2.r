#!/usr/bin/env Rscript

# Usage: 
# DECIPHER_DetectRepeat_MSA.2.r \
# --infile=~/EW_project/Kgroup/nt/K24316.fas.gz \
# --DetectRepeat=~/EW_project/TR02.16/DetectRepeats_nt/K24316.DetectRepeats.rds \
# --MSA=~/EW_project/TR02.16/MSA_nt/K24316.DetectRepeats.rds \
# --nt

# load libraries

library(optparse)

# specify options in a list
option_list = list(
	make_option(c("-i", "--infile"), type="character", default="NA", help="input seq (Required)", metavar="filename"),
	make_option(c("-r", "--DetectRepeat"), type="character", default="NA", help="output rds (Required)", metavar="filename"),
	make_option(c("-a", "--MSA"), type="character", default="NA", help="output rds (Required)", metavar="filename"),
	make_option(c("-n", "--nt"), action="store_true", default=FALSE, help="input is DNA seq [default=%default]")
); 

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(usage = "%prog [options] [--infile=FILENAME] [--DetectRepeat=FILENAME] [--MSA=FILENAME]", option_list=option_list)) 


######################## main ########################
library(DECIPHER, quietly = TRUE)

# read input seq file
if ( opt$nt ) {
	# print input filename
	message <- cat("readDNAStringSet", opt$infile)
	write(message, stderr())
    seq <- readDNAStringSet(opt$infile)
    write('AlignTranslation', stderr())
    result_MSA <- AlignTranslation(seq, processors=NULL, verbose=FALSE)
} else {
	# print input filename
	message <- cat("readAAStringSet", opt$infile)
	write(message, stderr())
    seq <- readAAStringSet(opt$infile)
    write('AlignSeqs', stderr())
    result_MSA <- AlignSeqs(seq, processors=NULL, verbose=FALSE)
}

# run DetectRepeats
result_DetectRepeats <- DetectRepeats(seq, processors=NULL, minScore = 0, allScores=TRUE, verbose=FALSE)
write('DetectRepeats', stderr())

# write output
saveRDS(result_MSA, file=opt$MSA, compress = TRUE)
saveRDS(result_DetectRepeats, file=opt$DetectRepeat, compress = TRUE)