require(dREG)

## Process command arguments
args <- commandArgs(trailingOnly=TRUE)

## Read PRO-seq data.
ps_plus_path  <- args[1]
ps_minus_path <- args[2]
outfile <- args[3];

cpu_cores <- as.integer(args[5])
if (is.na(cpu_cores)) cpu_cores <- 1;

gpu_id <- as.integer(args[6])
if (!is.na(gpu_id))
{
  library(Rgtsvm);
  ret <- selectGPUdevice( gpu_id );
}

use_rgtsvm=!is.na(gpu_id);

cat("------------ Parameters ------------- \n");
cat("Bigwig(plus):", ps_plus_path, "\n");
cat("Bigwig(minus):", ps_minus_path, "\n");
cat("Output:", outfile, "\n");
cat("dREG model:", args[4], "\n");
cat("CPU cores:", cpu_cores, "\n");
cat("GPU ID:", gpu_id, "\n");
cat("Using Rgtsvm:", use_rgtsvm, "\n");
cat("------------------------------------- \n");

if(!file.exists(ps_plus_path))
	stop( paste("Can't find the bigwig of plus strand(", ps_plus_path, ")"));
if(!file.exists(ps_minus_path))
	stop( paste("Can't find the bigwig of minus strand(", ps_minus_path, ")"));
if(!file.exists(args[4]))
	stop( paste("Can't find the SVR model(", args[4], ")"));


## Load the model.  Do this before loading ps_plus_path, just in case those are saved in the model file.
dreg_model <- args[4]
load(dreg_model) ## Should have (by default) gdm and asvm.


cat("[", as.character(Sys.time()), "] 1) Checking bigWig files.\n");
b1 <- check_bigwig(ps_plus_path, strand="+" );
b2 <- check_bigwig(ps_minus_path, strand="-" );
if( !b1 || !b2 )
{
    cat("Warning: bigWig files maybe not meet the requirements. See dREG requirement in https://github.com/Danko-Lab/dREG#data-preparation\n");
    stop("Stop");
}

cat("[", as.character(Sys.time()), "] 2) Scaning informative positions.\n");

## Now scan all positions in the genome ...
inf_positions <- get_informative_positions(ps_plus_path, ps_minus_path, depth= 0, step=50, use_ANDOR=TRUE, use_OR=FALSE) ## Get informative positions.

cat("Genome Loci=", NROW(inf_positions), "\n");

cat("[", as.character(Sys.time()), "] 3) Calling dREG prediction.\n");

t <- system.time( pred_val<- eval_reg_svm(gdm, asvm, inf_positions, ps_plus_path, ps_minus_path, batch_size= 50000, ncores=cpu_cores, use_rgtsvm = use_rgtsvm ) )

cat("[", as.character(Sys.time()), "] 4) Saving the result to compressed bed files.\n");

final_data <- data.frame(inf_positions, pred_val);
options("scipen"=100, "digits"=4)

write.table(final_data, file=paste(outfile, ".bed", sep=""), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

cat("Running time [User]:", t[1], "[System]:", t[2], "[Elapsed]:", t[3], "\n");
cat("Done!\n");


