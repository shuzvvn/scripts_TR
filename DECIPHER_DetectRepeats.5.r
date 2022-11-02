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
	make_option(c("-t", "--type"), type="character", default="tandem", help="Character string indicating the type of repeats to detect. This should be (an abbreviation of) one of 'tandem', 'interspersed', or 'both'. [default=%default]", metavar="string")
);

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(usage = "%prog [options] [--in_dataframe=FILENAME] [--row_index=NUMBER]", option_list=option_list)) 


######################## main ########################
suppressMessages(library(DECIPHER))
detectCores()
packageVersion("DECIPHER")

options(timeout=9999999) # default is 60 sec, will fail if the seq is too big

# read data frame that contains FTP address and genome names
genomes_df <- readRDS(opt$in_dataframe)
genomeIDs <- rownames(genomes_df)
genomeID_h <- genomeIDs[opt$row_index]
out_filename <- paste(genomeID_h, '.rds', sep='')
cat('\n', opt$row_index, genomeID_h, '\n')

fas.ftp <- genomes_df[genomeID_h, 'RefSeq.FTP']
if (fas.ftp == "") {
	fas.ftp <- genomes_df[genomeID_h, 'GenBank.FTP']
}
if (fas.ftp == "") {
	cat('\n', opt$row_index, genomeID_h, 'failed to get the FTP address.\n')
} else {
	# try downloading the cds sequence until succeed (at most 20 attempts)
	attempt <- 1
	fas.url <- paste(fas.ftp, "/", strsplit(fas.ftp, split = "/", fixed = TRUE)[[1]][10], "_protein.faa.gz", sep = "")
	cat("\nDownloading AA fasta file from: ", fas.url, '\n')
	attempt <- 1 # try downloading for at most 20 times
	while (!file.exists('protein.faa.gz') && attempt <= 20) {
		attempt <- attempt + 1
		try(
			download.file(fas.url, 'protein.faa.gz')
		)
	}
	if (!file.exists('protein.faa.gz')) {
		cat('\n', opt$row_index, genomeID_h, 'failed to download seq from FTP.\n')
	} else {
		cat('\nDONE downloading', opt$row_index, genomeID_h, 'from FTP.\n')
		end <- FALSE
		start_i <- 1
		while (end!=TRUE) {
			if (start_i==1) {
				# read in at most 100 seqs, starts from start_i (skip start_i-1)
				seq_h <- readAAStringSet('protein.faa.gz', nrec=100)
				cat('\nRunning DetectRepeats on CDS No.', start_i, 'to No.', start_i+length(seq_h)-1, '...\n\n')
				result_all <- DetectRepeats(seq_h, processors=8, verbose = TRUE, type = opt$type)
				start_i <- start_i + length(seq_h) # next round should starts from 101
			} else {
				# read in at most 100 seqs, starts from start_i
				seq_h <- readAAStringSet('protein.faa.gz', nrec=100, skip=start_i-1)
				if (length(seq_h)!=0) { # if there is seq
					cat('\nRunning DetectRepeats on CDS No.', start_i, 'to No.', start_i+length(seq_h)-1, '...\n\n')
					result_h <- DetectRepeats(seq_h, processors=8, verbose = TRUE, type = opt$type)
					result_h$Index <- result_h$Index + start_i-1
					result_all <- rbind(result_all, result_h)
					start_i <- start_i + length(seq_h) # next round should starts from
				} else {
					end <- TRUE
				}
			}
		}
		cat('\nDONE running DetectRepeats on', nrow(result_all), 'CDS...\n\n')
		# write output
		saveRDS(result_all, file=out_filename, compress = TRUE)
		cat('\nSave DetectRepeats result as', out_filename, '\n')
	}
}