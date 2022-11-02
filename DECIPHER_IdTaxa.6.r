# Usage: 
# Rscript DECIPHER_IdTaxa.5.r euk_1775.rds 1 KEGG_Eukaryotes_r95.RData
# load(url("http://www2.decipher.codes/Classification/TrainingSets/KEGG_Eukaryotes_r95.RData"))
# http://www2.decipher.codes/Classification/TrainingSets/KEGG_Animals_r95.RData

ARGS <- commandArgs(trailingOnly = TRUE)
in_dataframe <- ARGS[1] # in_dataframe
row_index <- as.integer(ARGS[2]) # row_index
TRS <- ARGS[3] # training set

print(ARGS)
print(in_dataframe)
print(row_index)

# load libraries
suppressMessages(library(DECIPHER))
detectCores()
packageVersion("DECIPHER")

options(timeout=9999999) # default is 60 sec, will fail if the genome is too big

######################## main ########################
# load training set
load(TRS, verbose = TRUE)

# load data
genomes_df <- readRDS(in_dataframe)
genomeIDs <- rownames(genomes_df)
i <- row_index
genomeID_h <- genomeIDs[i]
out_filename <- paste(genomeID_h, '.rds', sep='')
cat('\n', i, genomeID_h, '\n')
cat('\n', genomes_df[genomeID_h, 'Organism.Groups'], '\n')
cat("\nInput training set:", TRS, '\n')

fas.ftp <- genomes_df[genomeID_h, 'RefSeq.FTP']
if (fas.ftp == "") {
	fas.ftp <- genomes_df[genomeID_h, 'GenBank.FTP']
}
if (fas.ftp == "") {
	cat('\n', i, genomeID_h, 'failed to get the FTP address.\n')
} else {
	fas.url <- paste(fas.ftp, "/", strsplit(fas.ftp, split = "/", fixed = TRUE)[[1]][10], "_protein.faa.gz", sep = "")
	cat("\nDownloading AA fasta file from: ", fas.url, '\n')
	# try downloading the whole genome sequence until succeed
	attempt <- 1
	while (!file.exists('protein.faa.gz') && attempt <= 20) {
		attempt <- attempt + 1
		try(
			download.file(fas.url, 'protein.faa.gz')
		)
	}
	if (!file.exists('protein.faa.gz')) {
		cat('\n', i, genomeID_h, 'failed to download from FTP.\n')
	} else {
		cat('\nDONE downloading', i, genomeID_h, 'from FTP.\n')
		end <- FALSE
		start_i <- 1
		while (end!=TRUE) {
			if (start_i==1) {
				# read in at most 100 seqs, starts from start_i (skip start_i-1)
				seq_h <- readAAStringSet('protein.faa.gz', nrec=100)
				cat('\nRunning IdTaxa on CDS No.', start_i, 'to No.', start_i+length(seq_h)-1, '...\n\n')
				result_all <- IdTaxa(test = seq_h, trainingSet = trainingSet, threshold = 40, fullLength = 0.99, processors=8, verbose = TRUE)
				start_i <- start_i + length(seq_h) # next round should starts from 101
			} else {
				# read in at most 100 seqs, starts from start_i
				seq_h <- readAAStringSet('protein.faa.gz', nrec=100, skip=start_i-1)
				if (length(seq_h)!=0) { # if there is seq
					cat('\nRunning IdTaxa on CDS No.', start_i, 'to No.', start_i+length(seq_h)-1, '...\n\n')
					result_h <- IdTaxa(test = seq_h, trainingSet = trainingSet, threshold = 40, fullLength = 0.99, processors=8, verbose = TRUE)
					result_all <- c(result_all, result_h)
					start_i <- start_i + length(seq_h) # next round should starts from
				} else {
					end <- TRUE
				}
			}
		}
		cat('\nDONE running IdTaxa on', length(result_all), 'CDS...\n\n')
		# write output
		saveRDS(result_all, file=out_filename, compress = TRUE)
		cat('\nSave IdTaxa result as', out_filename, '\n')
	}
}