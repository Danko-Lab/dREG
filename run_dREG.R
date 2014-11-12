require(dREG)

## Process command arguments
args <- commandArgs(trailingOnly=TRUE)

## Read PRO-seq data.
ps_plus_path  <- args[1]
ps_minus_path <- args[2]
outfile <- args[3]
dreg_model <- args[4]
ncores <- as.integer(args[5])

load(dreg_model) ## Should have (by default) gdm and asvm.

## Now scan all positions in the genome ...
inf_positions <- get_informative_positions(ps_plus_path, ps_minus_path, depth= 0, step=50, use_ANDOR=TRUE, use_OR=FALSE) ## Get informative positions.

pred_val<- eval_reg_svm(gdm, asvm, inf_positions, ps_plus_path, ps_minus_path, batch_size= 10000, ncores=ncores)

final_data <- data.frame(inf_positions, pred_val)
options("scipen"=100, "digits"=4)
write.table(final_data, file=outfile, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")


