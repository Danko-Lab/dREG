require(dREG)
options("scipen"=100, "digits"=4)

## Process command arguments
args <- commandArgs(trailingOnly=TRUE)

## Load the dRGE model including two ojects 'asvm' and 'gdm'.
## Do this before loading ps_plus_path, just in case those are saved in the model file.
## Should have (by default) gdm and asvm.
load(args[4])

## Read arguments from thwe web page.
ps_plus_path  <- args[1]
ps_minus_path <- args[2]
## Read arguments from default parameters in run_dREG.sh
outfile <- args[3]

cpu_cores <- as.integer(args[5])
if (is.na(cpu_cores)) cpu_cores <- 1;

use_rgtsvm <- FALSE;
use_gpu <- toupper(as.character(args[6]))
if (!is.na(use_gpu) && (use_gpu=="GPU" || use_gpu=="TRUE") )
	use_rgtsvm <- TRUE;
if (!is.na(use_gpu) && (use_gpu=="FALSE") )
	use_rgtsvm <- FALSE;

gpu_cores <- as.integer(args[7])
if (is.na(gpu_cores)) gpu_cores <- 1;

cat("Bigwig(plus):", ps_plus_path, "\n");
cat("Bigwig(minus):", ps_minus_path, "\n");
cat("Output:", outfile, "\n");
cat("dREG model:", args[4], "\n");
cat("CPU cores:", cpu_cores, "\n");
cat("GPU:", use_rgtsvm, "\n");
cat("GPU cores:", gpu_cores, "\n");


if(!file.exists(ps_plus_path))
	stop( paste("Can't find the bigwig of plus strand(", ps_plus_path, ")"));
if(!file.exists(ps_minus_path))
	stop( paste("Can't find the bigwig of minus strand(", ps_minus_path, ")"));
if(!file.exists(args[4]))
	stop( paste("Can't find the SVR model(", args[4], ")"));

## Now scan all positions in the genome ...
cat("1) -------- Checking the informative positions\n");
load(args[4]);

cat("[", as.character(Sys.time()), "]", "Starting peak calling", "\n");
run.time <- system.time(r <- peak_calling( asvm, gdm, ps_plus_path, ps_minus_path, cpu_cores=cpu_cores, use_rgtsvm=use_rgtsvm, gpu_cores=gpu_cores));
cat("[", as.character(Sys.time()), "]", "Ending peak calling", "\n");

show( run.time/60 );

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

save(r, run.time, file=out.rdata );

#Rounding to 2 digits for predicted scores
r$infp_bed[,4] <- round( r$infp_bed[,4], 3 );
r$peak_bed[,4] <- round( r$peak_bed[,4], 3 );

make_index_gz( r$infp_bed, out.file1 );
make_index_gz( r$peak_bed, out.file2 );
make_index_gz( r$peak_bed[,c(1:4)], out.file3 );
make_index_gz( r$peak_broad, out.file4 );


#system( paste("tar -cvzf ", outfile, ".dREG.tar.gz", " ", outfile, ".dREG.*", sep="") );
#cat("Result:", paste(outfile, ".dREG.tar.gz", sep=""), "\n");

cat("Done!\n");
