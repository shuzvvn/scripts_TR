#!/usr/bin/Rscript

# Usage: 
# consist_DECIPHER.4.r \
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
	make_option(c("-k", "--Kgroups_uniq"), type="character", default="NA", help="Kgroups_uniq.rds (Required)", metavar="filename"),
    make_option(c("-r", "--DetectRepeat_dir"), type="character", default="NA", help="Directory path for Kgroup.DetectRepeat.rds (Required)", metavar="dirname"),
	make_option(c("-a", "--MSA_dir"), type="character", default="NA", help="Directory path for Kgroup.MSA.rds (Required)", metavar="dirname"),
	make_option(c("-s", "--score"), type="double", default=10, help="score cutoff (Required)", metavar="number"),
	make_option(c("-o", "--outfile"), type="character", default="NA", help="output tsv (Required)", metavar="filename"),
    make_option(c("-c", "--outcons"), type="character", default="NA", help="output cons (Required)", metavar="filename")
); 

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(usage = "%prog [options] [--DetectRepeat_dir=DIRNAME] [--Kgroups_uniq=FILENAME] [--MSA_dir=DIRNAME]  [--score=FLOAT] [--outfile=FILENAME] [--outcons=FILENAME]", option_list=option_list)) 


######################## main ########################

score_h <- opt$score

Kgroups_uniq <- readRDS(opt$Kgroups_uniq)


library(DECIPHER)
library(dplyr)
library(stringr)



total_n_seqs <- 0
total_n_seqs_r <- 0

# 3. repeats per seqs
total_n_r <- 0

# Consistency:
consistency <- c()

total_n_Kgroups <- 0 # number of Kgroups processed
total_n_Kgroups_r <- 0

for (Kgroup_h in Kgroups_uniq)
{
	total_n_Kgroups <- total_n_Kgroups + 1

    # get the filtered results
    DetectRepeat_filename <- paste(opt$DetectRepeat_dir, Kgroup_h, '.DetectRepeat.rds', sep='')
    DRresult_h <- readRDS(DetectRepeat_filename)
    DRresult_h <- DRresult_h[DRresult_h$Score >= opt$score,]

    # read MSA
    MSA_filename <- paste(opt$MSA_dir, Kgroup_h, '.MSA.rds', sep='')
    MSA_h <- readRDS(MSA_filename)
    MSA_len <- width(MSA_h)[1]
    n_seqs <- length(MSA_h) # seq in Kgroup
    n_seqs_r <- n_distinct(DRresult_h$Index) # seq have repeat 
    cat(paste(total_n_Kgroups, Kgroup_h, n_seqs, n_seqs_r, sep='\t'), '\n')

    total_n_seqs <- total_n_seqs + n_seqs
    total_n_seqs_r <- total_n_seqs_r + n_seqs_r
    
	if (nrow(DRresult_h)!=0) # ignore K group if no repeat detect
	{
		total_n_Kgroups_r <- total_n_Kgroups_r + 1
        total_n_r <- total_n_r + nrow(DRresult_h)
        # matrix of MSA, gap as NA
		mt_MSA <- matrix(, nrow = MSA_len, ncol = n_seqs)
		for (index_h in seq(n_seqs))
		{
			mt_MSA[,index_h] <- strsplit(toString(MSA_h[index_h]),"")[[1]]
		}
		mt_MSA[mt_MSA=='-'] <- NA

		# matrix of repeat loci in MSA
		mt_rep_loci <- matrix(0, nrow = MSA_len, ncol = n_seqs)
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
        # fit a function
        dataset <- cbind(pair_identity, pair_consistency)
        y.hat <- predict(lm(dataset$pair_consistency~dataset$pair_identity, data=dataset))
        corrected_pair_consistency <- pair_consistency / y.hat
        
        # report result
        cat(paste(pair_consistency, mean(corrected_pair_consistency), sep='\t'), '\n')
		consistency <- c(consistency, mean(corrected_pair_consistency))
        line <- paste(Kgroup_h, n_seqs, n_seqs_r, consistency, sep='\t')
        write(line, file=opt$outfile,append=TRUE)        
	}
}

cat(c(score_h, summary(consistency)[c(1:6)], total_n_seqs_r/total_n_seqs, total_n_Kgroups_r/total_n_Kgroups, total_n_r/total_n_seqs, '\n'), sep="\t", file=opt$outcons)