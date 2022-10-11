#!/usr/bin/env Rscript

# Usage: 
# count_cons_seqs.1.r --DR_filename=~/project/TR02/TR02.26/run1/aa_DetectRepeats/K24316.DetectRepeats.rds --MSA_filename=~/project/TR02/TR02.26/run1/aa_MSA/K24316.MSA.rds --score=15 --outfile=~/project/TR02/TR02.26/run1/aa_count_cons_seqs/K24316


# load libraries

library(optparse)

# specify options in a list
option_list = list(
	make_option(c("-r", "--DR_filename"), type="character", default="NA", help="File path for Kgroup.DetectRepeat.rds (Required)", metavar="filename"),
	make_option(c("-a", "--MSA_filename"), type="character", default="NA", help="File path for Kgroup.MSA.rds (Required)", metavar="filename"),
	make_option(c("-s", "--score"), type="double", default=15, help="score cutoff (Required)", metavar="number"),
	make_option(c("-o", "--outfile"), type="character", default="NA", help="output file (Required)", metavar="filename")
); 

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(usage = "%prog [options] [--DetectRepeat_dir=FILENAME] [--MSA_dir=FILENAME] [--score=FLOAT] [--outfile=FILENAME]", option_list=option_list)) 


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
get_mt_rep_loci <- function(DRresult_h, mt_MSA, MSA_len, n_seqs) {
	mt_rep_loci <- mt_MSA
	mt_rep_loci[!is.na(mt_rep_loci)] <- 0
	for (row_h in seq(nrow(DRresult_h)))
	{
		index_h <- DRresult_h[row_h, 'Index']
		begin_h <- DRresult_h[row_h, 'Begin']
		end_h <- DRresult_h[row_h, 'End']

		# convert location in ori_seq to location in alignment(with "-")
		# marked as:
		# 1: if location is in repeat region
		# 0: if location is not in repeat region
		# NA: if is gap
		pos_ori <- 0
		for (pos_align in seq(MSA_len))
		{
			if (!is.na(mt_rep_loci[pos_align,index_h]))
			{
				pos_ori <- pos_ori + 1
				if (between(pos_ori, begin_h, end_h))
				{
					mt_rep_loci[pos_align, index_h] <- 1
				}
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


suppressMessages(library(DECIPHER))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))

# cutoff 15
DRresult_h <- readRDS(opt$DR_filename)
DRresult_h <- DRresult_h[DRresult_h$Score >= opt$score,]

# read MSA
MSA_h <- readRDS(opt$MSA_filename)
MSA_len <- width(MSA_h)[1]

n_seqs <- length(MSA_h) # number of seq in Kgroup
n_seqs_r <- n_distinct(DRresult_h$Index) # seq have repeat

consider_n_pairs <- 0
has_cons_TR_n_pairs <- 0

# out_mt <- matrix(0, nrow = n_seqs, ncol = 3)
# colnames(out_mt) <- c('index', 'in_pairs', 'cons_repeat')
# out_mt[,'index'] <- c(1:n_seqs)

if (nrow(DRresult_h)!=0) # ignore K group if no repeat detect
{
	# MSA
	mt_MSA <- get_mt_MSA(MSA_h, MSA_len, n_seqs)
	# matrix of repeat loci in MSA
	mt_rep_loci <- get_mt_rep_loci(DRresult_h, mt_MSA, MSA_len, n_seqs)
	
	# pairwise consistency for positive locations ("1"), mask gaps (NA)
	pairs <- combn(seq(n_seqs), 2)
	for (pair_h in seq(ncol(pairs)))
	{
		# only consider pairs with PID between 50%-60%
		PID <- 1-DistanceMatrix(MSA_h[pairs[,pair_h]], type="dist", verbos=FALSE)[1]
		if (!is.na(PID))
		{
			if (PID >= 0.5 && PID <= 0.6)
			{
				#out_mt[pairs[,pair_h],'in_pairs'] <- 1 # seqs in pair
				pair_df <- as.data.frame(mt_rep_loci[,pairs[,pair_h]])
				if ( !all(is.na(pair_df[pair_df[,1] == 1, 2])) && !all(is.na(pair_df[pair_df[,2] == 1, 1]))) # exclude if repeat positions of one seq are all NA in another
				{
					# number of pair considered +=1
					consider_n_pairs <- consider_n_pairs + 1
					# remove rows if the alignment position is NA
					pair_df <- sapply(pair_df[rowSums(is.na(pair_df))==0,], as.numeric)
					if ( sum(pair_df[,1]*pair_df[,2]) > 0 ) # if any shared repeat position
					{
						has_cons_TR_n_pairs <- has_cons_TR_n_pairs + 1
					} 
				}
			}
		}
	}
	out_text <- c(ncol(pairs), consider_n_pairs, has_cons_TR_n_pairs)
} else {
	# MSA
	mt_MSA <- get_mt_MSA(MSA_h, MSA_len, n_seqs)
	pairs <- combn(seq(n_seqs), 2)
	for (pair_h in seq(ncol(pairs)))
	{
		# only consider pairs with PID between 50%-60%
		PID <- 1-DistanceMatrix(MSA_h[pairs[,pair_h]], type="dist", verbos=FALSE)[1]
		if (!is.na(PID))
		{
			if (PID >= 0.5 && PID <= 0.6)
			{
				consider_n_pairs <- consider_n_pairs + 1
			}
		}
	}
	out_text <- c(ncol(pairs), consider_n_pairs, has_cons_TR_n_pairs)
}

cat(out_text, file = opt$outfile, sep = "\t")