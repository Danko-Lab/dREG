require(dREG)
options("scipen"=100, "digits"=4)

## Process command arguments
args <- commandArgs(trailingOnly=TRUE);

## Read arguments from thwe web page.
ps_plus_path  <- args[1]
ps_minus_path <- args[2]
## Read arguments from default parameters in run_dREG.sh
outfile <- args[3]

cpu_cores <- as.integer(args[5])
if (is.na(cpu_cores)) cpu_cores <- 1;

use_rgtsvm <- FALSE;
gpu_id <- as.integer(args[6])
if (!is.na(gpu_id)) use_rgtsvm<-TRUE

cat("------------ Parameters ------------- \n");
cat("Bigwig(plus):", ps_plus_path, "\n");
cat("Bigwig(minus):", ps_minus_path, "\n");
cat("Output:", outfile, "\n");
cat("dREG model:", args[4], "\n");
cat("CPU cores:", cpu_cores, "\n");
cat("GPU ID:", gpu_id, "\n");
cat("Using Rgtsvm:", use_rgtsvm, "\n");
cat("-------------------------------------\n ");

if(!file.exists(ps_plus_path))
	stop( paste("Can't find the bigwig of plus strand(", ps_plus_path, ")"));
if(!file.exists(ps_minus_path))
	stop( paste("Can't find the bigwig of minus strand(", ps_minus_path, ")"));
if(!file.exists(args[4]))
	stop( paste("Can't find the SVR model(", args[4], ")"));

## Load the dRGE model including two ojects 'asvm' and 'gdm'.
## Do this before loading ps_plus_path, just in case those are saved in the model file.
## Should have (by default) gdm and asvm.
load(args[4]);

if(gpu_id>0)
{
	library(Rgtsvm);
	ret <- selectGPUdevice(gpu_id);
}

cat("[", as.character(Sys.time()), "] 1) Checking bigWig files.\n");
b1 <- check_bigwig(ps_plus_path, strand="+" );
b2 <- check_bigwig(ps_minus_path, strand="-" );
if( !b1 || !b2 )
{
    cat("Warning: bigWig files maybe not meet the requirements. See dREG requirement in https://github.com/Danko-Lab/dREG#data-preparation\n");
    stop("Stop");
}

## Now scan all positions in the genome ...
cat("[", as.character(Sys.time()), "] 2) Starting peak calling.\n");
run.time <- system.time(r <- peak_calling( asvm, gdm, ps_plus_path, ps_minus_path, cpu_cores=cpu_cores, use_rgtsvm=use_rgtsvm, gpu_cores=1));

r$run.time = run.time;

## Now scan all positions in the genome ...
cat("[", as.character(Sys.time()), "] 3) Saving the result to compressed bed files.\n");

out.file0 <- paste(outfile, "dREG.raw.peak.bed", sep=".")
out.file1 <- paste(outfile, "dREG.infp.bed", sep=".")
out.file2 <- paste(outfile, "dREG.peak.full.bed", sep=".")
out.file3 <- paste(outfile, "dREG.peak.score.bed", sep=".")
out.file4 <- paste(outfile, "dREG.peak.broad.bed", sep=".")
out.rdata <- paste(outfile, "dREG.rdata", sep=".")

make_index_gz<-function( df_bed, out_file)
{
	file.tmp <- tempfile(fileext=".bed");
	write.table( df_bed, file=file.tmp, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t");
	system(paste( "sort-bed ", file.tmp, " | bgzip -f > ", out_file, ".gz", sep="") );
	system(paste( "tabix -f -p bed ", out_file, ".gz", sep="") );
	unlink(file.tmp)
}

save(r, file=out.rdata );

#Rounding to 2 digits for predicted scores
r$infp_bed[,4] <- round( r$infp_bed[,4], 3 );
r$peak_bed[,4] <- round( r$peak_bed[,4], 3 );

make_index_gz( r$raw_peak[,-2], out.file0 );
make_index_gz( r$infp_bed, out.file1 );
make_index_gz( r$peak_bed, out.file2 );
make_index_gz( r$peak_bed[,c(1:4)], out.file3 );
make_index_gz( r$peak_broad, out.file4 );

#system( paste("tar -cvzf ", outfile, ".dREG.tar.gz", " ", outfile, ".dREG.*", sep="") );
#cat("Result:", paste(outfile, ".dREG.tar.gz", sep=""), "\n");

cat("Done!\n");