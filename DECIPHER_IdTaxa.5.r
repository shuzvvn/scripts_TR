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

options(timeout=3000000) # default is 60 sec, will fail if the genome is too big

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
cat("\nInput training set:", TRS)

seq <- NULL
fas.ftp <- genomes_df[genomeID_h, 'RefSeq.FTP']
if (fas.ftp == "") {
	fas.ftp <- genomes_df[genomeID_h, 'GenBank.FTP']
}
if (fas.ftp == "") {
	cat('\n', i, genomeID_h, 'failed to get the FTP address.\n')
} else {
	fas.url <- paste(fas.ftp, "/", strsplit(fas.ftp, split = "/", fixed = TRUE)[[1]][10], "_protein.faa.gz", sep = "")
    cat("\nreadAAStringSet from: ", fas.url, '\n')
	# try downloading the whole genome sequence until succeed
	attempt <- 1
	while (is.null(seq) && attempt <= 20) {
		attempt <- attempt + 1
		try(
			seq <- readAAStringSet(fas.url)
		)
	}
	if (is.null(seq)) {
		cat('\n', i, genomeID_h, 'failed to download from FTP.\n')
	} else {
        cat('\nRunning IdTaxa on', length(seq), ' CDS...\n\n')
        # run IdTaxa
        result_IdTaxa <- IdTaxa(test = seq, trainingSet = trainingSet, threshold = 40, fullLength = 0.99, processors = 8, verbose = TRUE)
        cat('\nIdTaxa DONE!', length(seq), '\n\n')
        # write output
        saveRDS(result_IdTaxa, file=out_filename, compress = TRUE)
	}
}