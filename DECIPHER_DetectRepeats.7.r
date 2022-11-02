# Usage: 
# Rscript DECIPHER_DetectRepeats.7.r euk_1775_fail.rds 1

ARGS <- commandArgs(trailingOnly = TRUE)
in_dataframe <- ARGS[1] # in_dataframe
row_index <- as.integer(ARGS[2]) # row_index

print(ARGS)
print(in_dataframe)
print(row_index)

# load libraries
suppressMessages(library(DECIPHER))
detectCores()
packageVersion("DECIPHER")

options(timeout=999999999) # default is 60 sec, will fail if the genome is too big

######################## main ########################
# read data frame that contains FTP address and genome names
genomes_df <- readRDS(in_dataframe)
genomeIDs <- rownames(genomes_df)
genomeID_h <- genomeIDs[row_index]
out_filename <- paste(genomeID_h, '.rds', sep='')
cat('\n', row_index, genomeID_h, '\n')

fas.ftp <- genomes_df[genomeID_h, 'RefSeq.FTP']
if (fas.ftp == "") {
	fas.ftp <- genomes_df[genomeID_h, 'GenBank.FTP']
}
if (fas.ftp == "") {
	cat('\n', row_index, genomeID_h, 'failed to get the FTP address.\n')
} else {
	# try downloading the cds sequence until succeed (at most 20 attempts)
	attempt <- 1
	fas.url <- paste(fas.ftp, "/", strsplit(fas.ftp, split = "/", fixed = TRUE)[[1]][10], "_protein.faa.gz", sep = "")
	cat("\nDownloading AA fasta file from: ", fas.url, '\n')
	attempt <- 1 # try downloading for at most 20 times
	while (!file.exists('protein.faa.gz') && attempt <= 5) {
		attempt <- attempt + 1
		try(
			download.file(fas.url, 'protein.faa.gz')
		)
	}
	if (!file.exists('protein.faa.gz')) {
		cat('\n', row_index, genomeID_h, 'failed to download seq from FTP.\n')
	} else {
		cat('\nDONE downloading', row_index, genomeID_h, 'from FTP.\n')
		end <- FALSE
		start_i <- 1
		while (end!=TRUE) {
			if (start_i==1) {
				# read in at most 100 seqs, starts from start_i (skip start_i-1)
				seq_h <- readAAStringSet('protein.faa.gz', nrec=100)
                seq_h <- RemoveGaps(seq_h, "all")
				cat('\nRunning DetectRepeats on CDS No.', start_i, 'to No.', start_i+length(seq_h)-1, '...\n\n')
				result_all <- DetectRepeats(seq_h, processors=8, verbose = TRUE, type = "tandem")
				start_i <- start_i + length(seq_h) # next round should starts from 101
			} else {
				# read in at most 100 seqs, starts from start_i
				seq_h <- readAAStringSet('protein.faa.gz', nrec=100, skip=start_i-1)
                seq_h <- RemoveGaps(seq_h, "all")
				if (length(seq_h)!=0) { # if there is seq
					cat('\nRunning DetectRepeats on CDS No.', start_i, 'to No.', start_i+length(seq_h)-1, '...\n\n')
					result_h <- DetectRepeats(seq_h, processors=8, verbose = TRUE, type = "tandem")
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