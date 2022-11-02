# Usage: 
# Rscript DECIPHER_IdTaxa_org.1.r pro4077_SSUs.rds GTDB_r207-mod_April2022.RData pro4077_SSUs_IdTaxa_org.rds

ARGS <- commandArgs(trailingOnly = TRUE)
in_SSUs <- ARGS[1] # in_SSUs
TRS <- ARGS[2] # training set
out_filename <- ARGS[3] # out filename


print(ARGS)
print(in_SSUs)

# load libraries
suppressMessages(library(DECIPHER))
detectCores()
packageVersion("DECIPHER")

######################## main ########################
# load training set
load(TRS, verbose = TRUE)

# load data
SSUs <- readRDS(in_SSUs)
cat('\nRead in', length(SSUs), 'seqs.\n')
cat("\nInput training set:", TRS, '\n')

cat('\nRunning IDTAXA Classify Organisms on input SSUs.\n')
ids <- IdTaxa(test = SSUs, trainingSet = trainingSet, threshold = 30, processors=8, verbose = TRUE)
cat('\nDONE running IdTaxa on', length(ids), 'SSUs...\n\n')

# write output
saveRDS(ids, file=out_filename, compress = TRUE)
cat('\nSave IdTaxa result as', out_filename, '\n')