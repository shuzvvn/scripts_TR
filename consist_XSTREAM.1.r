#!/usr/bin/env Rscript

# Usage: 
# consist_XSTREAM.1.r \
# --Kgroups_uniq=~/project/TR02/TR02.21/rds/Kgroups_uniq.rds \
# --XSTREAM_dir=~/project/TR02/TR02.22/aa_DetectRepeats/ \
# --MSA_dir=~/project/TR02/TR02.21/aa_MSA/ \
# --out_tsv=out.tsv \
# --out_mean_cons=out.mean.cons \

# load libraries

library(optparse)

# specify options in a list
option_list = list(
	make_option(c("-k", "--Kgroups_uniq"), type="character", default="NA", help="Kgroups_uniq.rds (Required)", metavar="filename"),
	make_option(c("-r", "--XSTREAM_dir"), type="character", default="NA", help="Directory path for Kgroup.xls (Required)", metavar="dirname"),
	make_option(c("-a", "--MSA_dir"), type="character", default="NA", help="Directory path for Kgroup.MSA.rds (Required)", metavar="dirname"),
	make_option(c("-s", "--score"), type="double", default=10, help="score cutoff (Required)", metavar="number"),
	make_option(c("-o", "--out_tsv"), type="character", default="NA", help="output tsv (Required)", metavar="filename"),
	make_option(c("-c", "--out_mean_cons"), type="character", default="NA", help="output cons (Required)", metavar="filename")
); 

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(usage = "%prog [options] [--Kgroups_uniq=FILENAME] [--XSTREAM_dir=DIRNAME] [--MSA_dir=DIRNAME] [--score=FLOAT] [--out_tsv=FILENAME] [--out_mean_cons=FILENAME]", option_list=option_list)) 


######################## functions ########################
# matrix of MSA, gap as NA
get_mt_MSA <- function(MSA_h, MSA_len, n_seqs) {
    # matrix of MSA, gap as NA
    mt_MSA <- matrix(, nrow = MSA_len, ncol = n_seqs)
    for (index_h in seq(n_seqs))
    {
        mt_MSA[,index_h] <- strsplit(toString(MSA_h[index_h]),"")[[1]]
    }
    mt_MSA[mt_MSA=='-'] <- NA
    return(mt_MSA)
}

# matrix of repeat loci in MSA
get_mt_rep_loci <- function(DRresult_h, MSA_h, mt_MSA, MSA_len, n_seqs) {
    mt_rep_loci <- matrix(0, nrow = MSA_len, ncol = n_seqs)
    for (row_h in seq(nrow(DRresult_h)))
    {
        index_h <- which(names(MSA_h) == DRresult_h[row_h, 'identifier'])
        begin_h <- DRresult_h[row_h, 'start']
        end_h <- DRresult_h[row_h, 'end']

        # convert location in ori_seq to location in alignment(with "-")
        # marked as:
        # 1: if location is in repeat region
        # 0: if location is not in repeat region
        # NA: if is gap
        pos_ori <- 0
        for (pos_align in seq(MSA_len))
        {
            if (!is.na(mt_MSA[pos_align,index_h]))
            {
                pos_ori <- pos_ori + 1
                if (between(pos_ori, begin_h, end_h))
                {
                    mt_rep_loci[pos_align, index_h] <- 1
                }
            } else {
                mt_rep_loci[pos_align, index_h] <- NA
            }
        }
    }
    return(mt_rep_loci)
}

# fraction of seq is in repeat region
get_frac_repeat_h <- function(mt_rep_loci) {
    n_repeat_loci <- colSums(mt_rep_loci==1,na.rm=TRUE)
    n_not_repeat_loci <- colSums(mt_rep_loci==0,na.rm=TRUE)
    return(n_repeat_loci / (n_repeat_loci+n_not_repeat_loci))
}

# pairwise consistency for positive locations ("1"), mask gaps (NA)
get_pair_consistency <- function(n_seqs, MSA_h, mt_rep_loci) {
    pairs <- combn(seq(n_seqs), 2)
    pair_consistency <- c()
    for (pair_h in seq(ncol(pairs)))
    {
        # only consider pairs with PID between 50%-60%
        PID <- 1-DistanceMatrix(MSA_h[pairs[,pair_h]], type="dist", verbos=FALSE)[1]
        #print(DistanceMatrix(MSA_h[pairs[,pair_h]], type="dist", verbos=FALSE))
        if (!is.na(PID))
        {
            if (PID >= 0.5 && PID <= 0.6)
            {
                pair_df <- as.data.frame(mt_rep_loci[,pairs[,pair_h]])
                if ( !all(is.na(pair_df[pair_df[,1] == 1, 2])) && !all(is.na(pair_df[pair_df[,2] == 1, 1]))) # exclude if repeat positions of one seq are all NA in another
                {
                    # remove rows if the alignment position is NA
                    pair_df <- pair_df[rowSums(is.na(pair_df))==0,]
                    # consistency is 0 if no repeat region found in one or both seqs of the pair after masking gaps
                    if ( 0 %in% colSums(pair_df) )
                    {
                        pair_consistency_h <- 0
                    } else {
                        # consistency is sum(both 1)/sum(at least one is 1)
                        pair_consistency_h <- sum(pair_df[,1]*pair_df[,2])/sum(rowSums(pair_df)!=0)
                    }
                    pair_consistency <- c(pair_consistency, pair_consistency_h)
                }
            }
        }
    }
    return(pair_consistency)
}



######################## main ########################

score_h <- opt$score
Kgroups_uniq <- readRDS(opt$Kgroups_uniq)

library(DECIPHER)
library(dplyr)
library(stringr)

# Y: Consistency:
consistency_mean <- c()
# X: mean(fraction of seq is in repeat region) %
frac_repeat <- c()

line <- paste('KgroupID', 'n_seqs', 'n_seqs_r', 'n_pairs', 'mean_cons', 'weighted_mean_cons', sep='\t')
write(line, file=opt$out_tsv, append=FALSE)

for (Kgroup_h in Kgroups_uniq)
{
	# get the filtered results
	XSTREAM_filename <- paste(opt$XSTREAM_dir, Kgroup_h, '.xls', sep='')
    DRresult_h <- read.table(file=XSTREAM_filename, header=TRUE, sep = "\t")

	# read MSA
	MSA_filename <- paste(opt$MSA_dir, Kgroup_h, '.MSA.rds', sep='')
	MSA_h <- readRDS(MSA_filename)
	MSA_len <- width(MSA_h)[1]
	
	n_seqs <- length(MSA_h) # number of seq in Kgroup
	n_seqs_r <- n_distinct(DRresult_h$identifier) # seq have repeat

	if (nrow(DRresult_h)!=0) # ignore K group if no repeat detect
	{
        mt_MSA <- get_mt_MSA(MSA_h, MSA_len, n_seqs)
        mt_rep_loci <- get_mt_rep_loci(DRresult_h, MSA_h, mt_MSA, MSA_len, n_seqs) 
        frac_repeat_h <- get_frac_repeat_h(mt_rep_loci)
        frac_repeat <- c(frac_repeat, frac_repeat_h)
        pair_consistency <- get_pair_consistency(n_seqs, MSA_h, mt_rep_loci)

		# report result
		# Weight the mean of a Kgroup base on # of pairs be considered
		if ( length(pair_consistency) > 0 )
		{
			mean_h <- mean(pair_consistency)
			weighted_mean_h <- mean_h*length(pair_consistency)
			consistency_mean <- c(consistency_mean, weighted_mean_h)
			line <- paste(Kgroup_h, n_seqs, n_seqs_r, length(pair_consistency), mean_h, weighted_mean_h, sep='\t')
			write(line, file=opt$out_tsv, append=TRUE)
		} else {
			line <- paste(Kgroup_h, n_seqs, n_seqs_r, 0, NA, NA, sep='\t')
			write(line, file=opt$out_tsv, append=TRUE)
		}
	} else {
		# if no repeat detect, add zeros to frac_repeat
		frac_repeat <- c(frac_repeat, integer(n_seqs))
		line <- paste(Kgroup_h, n_seqs, n_seqs_r, NA, NA, NA, sep='\t')
		write(line, file=opt$out_tsv, append=TRUE)
	}
}

# X-axis: mean(fraction of seq is in repeat region) %
mean_frac_repeat <- mean(frac_repeat)
if ( length(consistency_mean) > 0 )
{
	cat(c(score_h, summary(consistency_mean)[c(1:6)], mean_frac_repeat, '\n'), sep="\t", file=opt$out_mean_cons)
} else {
	cat(c(score_h, integer(6), mean_frac_repeat, '\n'), sep="\t", file=opt$out_mean_cons)
}
