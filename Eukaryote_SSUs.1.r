#!/usr/bin/env Rscript
# load libraries
library(optparse)

# specify options in a list
option_list = list(
	make_option(c("-n", "--start_index"), type="double", default=1, help="row number to start with", metavar="number"),
	make_option(c("-o", "--out"), type="character", default="NA", help="output .rds (Required)", metavar="filename")
); 

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(usage = "%prog [options] [--start_index=NUMBER] [--out=FILENAME]", option_list=option_list)) 


library(DECIPHER)

# Eukaryotes
euk_1775 <- readRDS('euk_1775.rds')
genomeIDs <- rownames(euk_1775)
data(NonCodingRNA_Eukarya)
x <- NonCodingRNA_Eukarya

Euk_failed <- c()
Euk_SSUs <- NULL
n <- opt$start_index

cat('Batch row number', n, 'to', n+24, '\n')

for (i in c(n:(n+24)))
{
	genomeID_h <- genomeIDs[i]
	cat(i, genomeID_h, '\n')
	fas.ftp <- euk_1775[genomeID_h, 'RefSeq.FTP']
	if (fas.ftp == "") {
		fas.ftp <- euk_1775[genomeID_h, 'GenBank.FTP']
	}
	if (fas.ftp == "") {
		Euk_failed <- c(Euk_failed, genomeID_h)
	} else {
		fas.url <- paste(fas.ftp, "/", strsplit(fas.ftp, split = "/", fixed = TRUE)[[1]][10], "_genomic.fna.gz", sep = "")
		# try downloading the whole genome sequence until succeed
		genome <- NULL
		attempt <- 1
		while (is.null(genome) && attempt <= 5) {
			attempt <- attempt + 1
			try(
				genome <- readDNAStringSet(fas.url)
			)
		}
		if (is.null(genome)) {
			Euk_failed <- c(Euk_failed, genomeID_h)
		} else {
			cat('n_chr:', length(genome), '\n')
			best_score <- 0
			for (i_chr in seq_along(genome)) # loop through chromosomes, find the SSU with the highest score
			{
				#chr_h <- genome[i_chr]
				print(names(genome[i_chr]))
				z <- FindNonCoding(x['rRNA_18S-RF01960'], genome[i_chr], processors=1) # find 16S seq
				print(z)
				if ( length(z) > 0 && max(z[, 'TotalScore']) > best_score) {
					genes <- ExtractGenes(z[which.max(z[, 'TotalScore'])], genome[i_chr], type="RNAStringSet") # get the ones with highest score
					best_score <- max(z[, 'TotalScore'])
				} 
			}
			if (best_score==0) {
				Euk_failed <- c(Euk_failed, genomeID_h)
			} else {
				SSU_h <- RNAStringSet(unname(as.character(genes[1]))) # pick the first one if there are multiple with the same score
				names(SSU_h) <- genomeID_h
				if (is.null(Euk_SSUs)) {
					Euk_SSUs <- SSU_h
				} else {
					Euk_SSUs <- c(Euk_SSUs, SSU_h)
				}
			}
		}

	}
}

cat('Failed Genomes:', Euk_failed, '\n')
saveRDS(Euk_SSUs, file = opt$out, compress = TRUE)
