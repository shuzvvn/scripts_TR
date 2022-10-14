#!/usr/bin/env Rscript

# Usage: 
# DECIPHER_DetectRepeats.4.r \
# --in_dataframe=euk_1775_fail.rds \
# --row_index=1 \
# --type="tandem"
# --nt

# load libraries
library(optparse)

# specify options in a list
option_list = list(
	make_option(c("-i", "--in_dataframe"), type="character", default="NA", help="input data frame (Required)", metavar="filename"),
	make_option(c("-r", "--row_index"), type="double", default=1, help="row number", metavar="number"),
	make_option(c("-n", "--nt"), action="store_true", default=FALSE, help="input is DNA seq [default=%default]"),
	make_option(c("-t", "--type"), type="character", default="tandem", help="Character string indicating the type of repeats to detect. This should be (an abbreviation of) one of 'tandem', 'interspersed', or 'both'. [default=%default]", metavar="string")
);

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(usage = "%prog [options] [--in_dataframe=FILENAME] [--row_index=NUMBER]", option_list=option_list)) 


######################## main ########################
suppressMessages(library(DECIPHER))
detectCores()
packageVersion("DECIPHER")

options(timeout=3000000) # default is 60 sec, will fail if the seq is too big

# read data frame that contains FTP address and genome names
genomes_df <- readRDS(opt$in_dataframe)
genomeIDs <- rownames(genomes_df)

genomeID_h <- genomeIDs[opt$row_index]

out_filename <- paste(genomeID_h, '.rds', sep='')

cat('\n', opt$row_index, genomeID_h, '\n')

fas.ftp <- genomes_df[genomeID_h, 'RefSeq.FTP']
seq <- NULL
if (fas.ftp == "") {
	fas.ftp <- genomes_df[genomeID_h, 'GenBank.FTP']
}
if (fas.ftp == "") {
	cat('\n', opt$row_index, genomeID_h, 'failed to get the FTP address.\n')
} else {
	# try downloading the cds sequence until succeed (at most 20 attempts)
	attempt <- 1
	if ( opt$nt ) {
		# print input filename
		fas.url <- paste(fas.ftp, "/", strsplit(fas.ftp, split = "/", fixed = TRUE)[[1]][10], "_cds_from_genomic.fna.gz", sep = "")
		cat("\nreadDNAStringSet from: ", fas.url, '\n')
		while (is.null(seq) && attempt <= 20) {
			attempt <- attempt + 1
			try(
				seq <- readDNAStringSet(fas.url)
			)
		}
	} else {
		# print input filename
		fas.url <- paste(fas.ftp, "/", strsplit(fas.ftp, split = "/", fixed = TRUE)[[1]][10], "_protein.faa.gz", sep = "")
		cat("\nreadAAStringSet from: ", fas.url, '\n')
		while (is.null(seq) && attempt <= 20) {
			attempt <- attempt + 1
			try(
				seq <- readAAStringSet(fas.url)
			)
		}
	}
	# Run DetectRepeats if seq download successfully
	if (is.null(seq)) {
		cat('\n', opt$row_index, genomeID_h, 'failed to download seq from FTP.\n')
	} else {
		cat('\nRunning DetectRepeats on', length(seq), ' CDS...\n\n')
		# run DetectRepeats
		result_DetectRepeats <- DetectRepeats(seq, processors=NULL, verbose=TRUE, type = opt$type)
		cat('\nDetectRepeats DONE!', length(seq), '\n\n')
		# write output
		saveRDS(result_DetectRepeats, file=out_filename, compress = TRUE)
	}
}