#!/usr/bin/Rscript

# Usage: 
# consist_DECIPHER.R \
# --list_DRresult=list_DRresult.rds \
# --Kgroups_uniq_subset=Kgroups_uniq_subset.rds \
# --list_MSA=list_MSA.rds \
# --sum_n_seqs=10000 \
# --n_subset=100 \
# --score=10 \
# --outfile=out.tsv

# load libraries

library(optparse)

# specify options in a list
option_list = list(
	make_option(c("-i", "--list_DRresult"), type="character", default="NA", help="list_DRresult.rds (Required)", metavar="filename"),
    make_option(c("-k", "--Kgroups_uniq_subset"), type="character", default="NA", help="Kgroups_uniq_subset.rds (Required)", metavar="filename"),
    make_option(c("-a", "--list_MSA"), type="character", default="NA", help="list_MSA.rds (Required)", metavar="filename"),
    make_option(c("-n", "--sum_n_seqs"), type="integer", default=10000, help="total number of sequences (Required)", metavar="number"),
    make_option(c("-m", "--n_subset"), type="integer", default=100, help="total number of K groups (Required)", metavar="number"),
    make_option(c("-s", "--score"), type="double", default=10, help="score cutoff (Required)", metavar="number"),
	make_option(c("-o", "--outfile"), type="character", default="NA", help="output tsv (Required)", metavar="filename")
); 

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(usage = "%prog [options] [--list_DRresult=FILENAME] [--Kgroups_uniq_subset=FILENAME] [--list_MSA=FILENAME] [--sum_n_seqs=INTERGER] [--n_subset=INTERGER] [--score=FLOAT] [--outfile=FILENAME]", option_list=option_list)) 


######################## main ########################

score_h <- opt$score
sum_n_seqs <- opt$sum_n_seqs
n_subset <- opt$n_subset
list_DRresult <- readRDS(opt$list_DRresult)
Kgroups_uniq_subset <- readRDS(opt$Kgroups_uniq_subset)
list_MSA <- readRDS(opt$list_MSA)

library(DECIPHER)
library(dplyr)

pairs <- combn(seq(n_subset), 2)

# get the filtered results
list_DRresult_h <- lapply(list_DRresult, function(x) filter(x, Score >= score_h))

# 1. number of seqs have repeat
frac_seqs_h <- sum(sapply(list_DRresult_h, function(x) n_distinct(x$Index))) / sum_n_seqs

# 2. number of Kgroups have repeat
frac_DB_h <- sum(sapply(list_DRresult_h, nrow) != 0) / n_subset

# 3. repeats per seqs
repeats_per_seqs_h <- sum(sapply(list_DRresult_h, nrow)) / sum_n_seqs
# Consistency: average identity for positives (%)
consistency <- c()

for (Kgroup_h in Kgroups_uniq_subset)
{
    DRresult_h <- list_DRresult_h[[Kgroup_h]]
    if (nrow(DRresult_h)!=0) # ignore K group if no repeat detect
    {
        # MSA
        MSA_h <- list_MSA[[Kgroup_h]]
        MSA_len <- width(MSA_h)[1]
        mt_rep_loc <- matrix(0, nrow = MSA_len, ncol = n_subset)
        for (row_h in seq(nrow(DRresult_h)))
        {
            index_h <- DRresult_h[row_h, 'Index']
            begin_h <- DRresult_h[row_h, 'Begin']
            end_h <- DRresult_h[row_h, 'End']
            message <- paste(c(score_h, Kgroup_h, row_h, index_h, begin_h, end_h), sep='\t')
            write(message, stderr())
            align_seq_h <- unlist(strsplit(toString(MSA_h[index_h]),""))
            # convert position in ori_seq to position in alignment(with "-")
            # marked as "1" if position in repeat region and not "-"
            pos_ori <- 0
            for (pos_align in seq(MSA_len))
            {
                if (align_seq_h[pos_align] != '-')
                {
                    pos_ori <- pos_ori + 1
                    if (between(pos_ori, begin_h, end_h))
                    {
                        mt_rep_loc[pos_align, index_h] <- 1
                    }
                }
            }
        }
        # pairwise identity for positive positions, ignore gaps
        idents <- c()
        for (pair_h in seq(ncol(pairs)))
        {
            index_1 <- pairs[1,pair_h]
            index_2 <- pairs[2,pair_h]
            n_rep_pos_1 <- sum(mt_rep_loc[,index_1])
            n_rep_pos_2 <- sum(mt_rep_loc[,index_2])
            n_rep_pos_both <- sum(mt_rep_loc[,index_1] & mt_rep_loc[,index_1]==mt_rep_loc[,index_2])
            # ident = 0 if no repeat region found
            if (n_rep_pos_1 == 0)
            {
                ident_1 <- 0
            } else {
                ident_1 <- n_rep_pos_both/n_rep_pos_1
            }
            if (n_rep_pos_2 == 0)
            {
                ident_2 <- 0
            } else {
                ident_2 <- n_rep_pos_both/n_rep_pos_2
            }
            
            idents <- c(idents, mean(c(ident_1, ident_2)))
        }
        consistency <- c(consistency, mean(idents))
    }
}
ave_cons_h <- mean(consistency)

cat(c(score_h, ave_cons_h, frac_seqs_h, frac_DB_h, repeats_per_seqs_h, '\n'), sep="\t", file=opt$outfile)