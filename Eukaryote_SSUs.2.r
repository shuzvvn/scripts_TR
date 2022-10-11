#!/usr/bin/env Rscript

# load libraries
library(optparse)

# specify options in a list
option_list = list(
    make_option(c("-i", "--in_dataframe"), type="character", default="NA", help="input data frame (Required)", metavar="filename"),
	make_option(c("-r", "--row_index"), type="double", default=1, help="row number", metavar="number")
); 

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(usage = "%prog [options] [--in_dataframe=FILENAME] [--start_index=NUMBER] [--out=FILENAME]", option_list=option_list)) 


library(DECIPHER)

# Eukaryotes
genomes_df <- readRDS(opt$in_dataframe)
genomeIDs <- rownames(genomes_df)
data(NonCodingRNA_Eukarya)
x <- NonCodingRNA_Eukarya

i <- opt$row_index

genomeID_h <- genomeIDs[i]
out_filename <- paste(genomeID_h, '_ssu.fasta', sep='')
cat(i, genomeID_h, '\n')
fas.ftp <- genomes_df[genomeID_h, 'RefSeq.FTP']
if (fas.ftp == "") {
    fas.ftp <- genomes_df[genomeID_h, 'GenBank.FTP']
}
if (fas.ftp == "") {
    cat(i, genomeID_h, 'FTP address not exist.\n')
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
        cat(i, genomeID_h, 'failed to download from FTP.\n')
    } else {
        cat('n_chr:', length(genome), '\n\n')
        best_score <- 0
        for (i_chr in seq_along(genome)) # loop through chromosomes, find the SSU with the highest score
        {
            print(names(genome[i_chr]))
            z <- FindNonCoding(x['rRNA_18S-RF01960'], genome[i_chr], processors=1) # find 16S seq
            print(z)
            if ( length(z) > 0 && max(z[, 'TotalScore']) > best_score) {
                genes <- ExtractGenes(z[which.max(z[, 'TotalScore'])], genome[i_chr], type="RNAStringSet") # get the ones with highest score
                best_score <- max(z[, 'TotalScore'])
            } 
        }
        if (best_score==0) {
            cat(i, genomeID_h, 'failed to find any 16S seq.\n')
        } else {
            SSU_h <- RNAStringSet(unname(as.character(genes[1]))) # pick the first one if there are multiple with the same score
            names(SSU_h) <- genomeID_h
            writeXStringSet(SSU_h, out_filename)
        }
    }
}
