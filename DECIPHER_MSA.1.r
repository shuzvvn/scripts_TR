#!/usr/bin/Rscript

# Usage: 
# DECIPHER_MSA.1.r \
# --infile=/mnt/c/Users/stc/project/TR02/source_data/aa.fasta/1A12A.fasta \
# --out_AlignSeqs=/mnt/c/Users/stc/project/TR02/TR02.09/1A12A.AlignSeqs.rds \
# --out_AlignTranslation=/mnt/c/Users/stc/project/TR02/TR02.09/1A12A.AlignTranslation.rds \
# --nt


# load libraries

library(optparse)

# specify options in a list
option_list = list(
    make_option(c("-i", "--infile"), type="character", default="NA", help="input seq (Required)", metavar="filename"),
    make_option(c("-s", "--out_AlignSeqs"), type="character", default="NA", help="output .rds (Required)", metavar="filename"),
    make_option(c("-t", "--out_AlignTranslation"), type="character", default="NA", help="output .rds", metavar="filename"),
    make_option(c("-n", "--nt"), action="store_true", default=FALSE, help="input is DNA seq [default=%default]")
); 

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(usage = "%prog [options] [--infile=FILENAME] [--out_AlignSeqs=FILENAME] [--out_AlignTranslation=FILENAME]", option_list=option_list)) 


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

# run AlignSeqs
result_AlignSeqs <- AlignSeqs(seq, processors=NULL)
result_AlignTranslation <- AlignTranslation(seq, processors=NULL)


# save results as RDS files
saveRDS(result_AlignSeqs, file=opt$out_AlignSeqs, compress = TRUE)
saveRDS(result_AlignTranslation, file=opt$out_AlignTranslation, compress = TRUE)