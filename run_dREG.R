require(dREG)

## Process command arguments
args <- commandArgs(trailingOnly=TRUE)

## Load the model.  Do this before loading ps_plus_path, just in case those are saved in the model file.
dreg_model <- args[4]
load(dreg_model) ## Should have (by default) gdm and asvm.

## Read PRO-seq data.
ps_plus_path  <- args[1]
ps_minus_path <- args[2]
outfile <- args[3]
ncores <- as.integer(args[5])
if (is.na(ncores)) ncores <- 1;

use_rgtsvm <- FALSE;
use_gpu <- toupper(as.character(args[6]))
if (!is.na(use_gpu) && ( use_gpu=="GPU" || use_gpu=="TRUE")  ) use_rgtsvm <- TRUE;

## Now scan all positions in the genome ...
inf_positions <- get_informative_positions(ps_plus_path, ps_minus_path, depth= 0, step=50, use_ANDOR=TRUE, use_OR=FALSE) ## Get informative positions.

cat("Genome Loci=", NROW(inf_positions), "\n");

t <- system.time( pred_val<- eval_reg_svm(gdm, asvm, inf_positions, ps_plus_path, ps_minus_path, batch_size= 50000, ncores=ncores, use_rgtsvm=use_rgtsvm) )

cat("Running time [User]:", t[1], "[System]:", t[2], "[Elapsed]:", t[3], "\n");

final_data <- data.frame(inf_positions, pred_val)
options("scipen"=100, "digits"=4)
write.table(final_data, file=outfile, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")


