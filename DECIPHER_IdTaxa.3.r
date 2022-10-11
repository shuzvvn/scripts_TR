#!/usr/bin/env Rscript

# Usage: 
# DECIPHER_IdTaxa.3.r /Users/stc/script/DECIPHER_IdTaxa.3.r --query=/Users/stc/project/TR02/source_data/GbBac_3822/faa/GCF_000007725.1_ASM772v1_protein.faa.gz --TRS=/Users/stc/project/TR02/TR02.41/KEGG_Gammaproteobacteria_r95.RData --DR=/Users/stc/project/TR02/TR02.33/aa_DetectRepeats/GCF_000007725.1_ASM772v1.DetectRepeats.rds --out=/Users/stc/project/TR02/TR02.42/GCF_000007725.IdTaxa.rds

# load libraries

library(optparse)

# specify options in a list
option_list = list(
	make_option(c("-i", "--query"), type="character", default="NA", help="input seq (Required)", metavar="filename"),
    make_option(c("-r", "--DR"), type="character", default="NA", help="input DetectRepeats result (.rds) (Required)", metavar="filename"),
	make_option(c("-t", "--TRS"), type="character", default="NA", help="training set .RData (Required)"),
    make_option(c("-o", "--out"), type="character", default="NA", help="output .rds (Required)", metavar="filename"),
	make_option(c("-s", "--threshold"), type="double", default=40, help="threshold", metavar="number"),
    make_option(c("-c", "--TR_cutoff"), type="double", default=30, help="threshold", metavar="number"),
    make_option(c("-n", "--nt"), action="store_true", default=FALSE, help="input is DNA seq [default=%default]")
); 

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(usage = "%prog [options] [--query=FILENAME] [--DR=FILENAME] [--TRS=FILENAME] [--out=FILENAME]", option_list=option_list)) 


######################## main ########################
suppressMessages(library(DECIPHER))

# source('utils.R')
# source('ScoreAlignment.R')

# load training set
load(opt$TRS, verbose = TRUE)
message <- cat("\nInput training set:", opt$TRS)
write(message, stderr())

# read input seq file
if ( opt$nt ) {
	# print input filename
	message <- cat("\nreadDNAStringSet", opt$query)
	write(message, stderr())
	seq <- translate(readDNAStringSet(opt$query), if.fuzzy="solve")
} else {
	# print input filename
	message <- cat("\nreadAAStringSet", opt$query)
	write(message, stderr())
	seq <- readAAStringSet(opt$query)
}

# read DetectRepeats results
DRresult_h <- readRDS(opt$DR)
DRresult_h <- DRresult_h[DRresult_h$Score >= opt$TR_cutoff,]
Index_h <- unique(DRresult_h$Index)
seq_TR <- seq[Index_h]

# run IdTaxa
write("\nIdTaxa\n", stderr())
result_IdTaxa <- IdTaxa(test = seq_TR, trainingSet = trainingSet, threshold = opt$threshold, processors = NULL, verbose = TRUE)

# write output
saveRDS(result_IdTaxa, file=opt$out, compress = TRUE)
