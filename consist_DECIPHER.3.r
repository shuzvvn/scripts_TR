#!/usr/bin/Rscript

# Usage: 
# consist_DECIPHER.3.r \
# --list_DRresult=list_DRresult.rds \
# --Kgroups_uniq=Kgroups_uniq.rds \
# --list_MSA=list_MSA.rds \
# --sum_n_seqs=10000 \
# --n_subset=100 \
# --score=10 \
# --outfile=out.tsv \
# --consistency=out.cons

# load libraries

library(optparse)

# specify options in a list
option_list = list(
	make_option(c("-i", "--list_DRresult"), type="character", default="NA", help="list_DRresult.rds (Required)", metavar="filename"),
	make_option(c("-k", "--Kgroups_uniq"), type="character", default="NA", help="Kgroups_uniq.rds (Required)", metavar="filename"),
	make_option(c("-a", "--list_MSA"), type="character", default="NA", help="list_MSA.rds (Required)", metavar="filename"),
	make_option(c("-n", "--sum_n_seqs"), type="integer", default=10000, help="total number of sequences (Required)", metavar="number"),
	make_option(c("-m", "--n_subset"), type="integer", default=100, help="total number of K groups (Required)", metavar="number"),
	make_option(c("-s", "--score"), type="double", default=10, help="score cutoff (Required)", metavar="number"),
	make_option(c("-o", "--outfile"), type="character", default="NA", help="output tsv (Required)", metavar="filename"),
	make_option(c("-c", "--consistency"), type="character", default="NA", help="output consistency (Required)", metavar="filename")
); 

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(usage = "%prog [options] [--list_DRresult=FILENAME] [--Kgroups_uniq=FILENAME] [--list_MSA=FILENAME] [--sum_n_seqs=INTERGER] [--n_subset=INTERGER] [--score=FLOAT] [--outfile=FILENAME] [--consistency=FILENAME]", option_list=option_list)) 


######################## main ########################

score_h <- opt$score
sum_n_seqs <- opt$sum_n_seqs
n_subset <- opt$n_subset
list_DRresult <- readRDS(opt$list_DRresult)
Kgroups_uniq <- readRDS(opt$Kgroups_uniq)
list_MSA <- readRDS(opt$list_MSA)

library(DECIPHER)
library(dplyr)
library(stringr)

# get the filtered results
list_DRresult_h <- lapply(list_DRresult, function(x) filter(x, Score >= score_h))

# 1. number of seqs have repeat
frac_seqs_h <- sum(sapply(list_DRresult_h, function(x) n_distinct(x$Index))) / sum_n_seqs

# 2. number of Kgroups have repeat
frac_DB_h <- sum(sapply(list_DRresult_h, nrow) != 0) / n_subset

# 3. repeats per seqs
repeats_per_seqs_h <- sum(sapply(list_DRresult_h, nrow)) / sum_n_seqs

# Consistency:
consistency <- c()

n <- 0 # number of Kgroups processed
for (Kgroup_h in Kgroups_uniq)
{
	n <- n + 1
	cat(paste(n, Kgroup_h, sep='\t'), '\n')
	DRresult_h <- list_DRresult_h[[Kgroup_h]]
	if (nrow(DRresult_h)!=0) # ignore K group if no repeat detect
	{
		# read MSA
		MSA_h <- list_MSA[[Kgroup_h]]
		MSA_len <- width(MSA_h)[1]
		n_seqs <- length(MSA_h)
		
		# matrix of MSA, gap as NA
		mt_MSA <- matrix(, nrow = MSA_len, ncol = n_seqs)
		for (index_h in seq(n_seqs))
		{
			mt_MSA[,index_h] <- strsplit(toString(MSA_h[index_h]),"")[[1]]
		}
		mt_MSA[mt_MSA=='-'] <- NA
		print('MSA matrix done')

		# matrix of repeat loci in MSA
		mt_rep_loci <- matrix(0, nrow = MSA_len, ncol = n_seqs)
		for (row_h in seq(nrow(DRresult_h)))
		{
			index_h <- DRresult_h[row_h, 'Index']
			begin_h <- DRresult_h[row_h, 'Begin']
			end_h <- DRresult_h[row_h, 'End']
			#align_seq_h <- strsplit(toString(MSA_h[index_h]),"")[[1]]
			#align_seq_h <- mt_MSA[,index_h]

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
		print('Repeat matrix done')

		# pairwise consistency for positive locations ("1"), mask gaps (NA)
		pairs <- combn(seq(n_seqs), 2)
		pair_consistency <- c()
		pair_identity <- c()
		for (pair_h in seq(ncol(pairs)))
		{
            pair_df <- as.data.frame(mt_rep_loci[,pairs[,pair_h]])
            if ( !all(is.na(pair_df[pair_df[,1] == 1, 2])) && !all(is.na(pair_df[pair_df[,2] == 1, 1]))) # exclude if repeat positions of one seq are all NA in another
            {
                # calculate PID2: 100 * (identical positions) / (aligned positions)
	            identical_positions <- sum(mt_MSA[,pairs[1,pair_h]] == mt_MSA[,pairs[2,pair_h]], na.rm=T)
	            aligned_positions <- (as.integer(!is.na(mt_MSA[,pairs[1,pair_h]])) %*% as.integer(!is.na(mt_MSA[,pairs[2,pair_h]])))[1] # exclude aligned gaps
	            PID_h = identical_positions/aligned_positions

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
                pair_identity <- c(pair_identity, PID_h)
            }
		}
		#cat(c(Kgroup_h, mean(pair_consistency), '\n'), sep='\t')
		consistency <- c(consistency, mean(pair_consistency))
		write.table(cbind(pair_identity, pair_consistency), file=paste("plot", Kgroup_h, score_h, "tsv", sep = "."), row.names=FALSE, sep="\t")
	}
}
summary_cons_h <- summary(consistency)[c(1:6)]

cat(c(score_h, summary_cons_h, frac_seqs_h, frac_DB_h, repeats_per_seqs_h, '\n'), sep="\t", file=opt$outfile)

write(paste(consistency, sep="\n"), file=opt$consistency)