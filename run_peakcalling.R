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

ncores <- as.integer(args[5])
if (is.na(ncores)) ncores <- 1;

use_rgtsvm <- FALSE;
use_gpu <- toupper(as.character(args[6]))
if (!is.na(use_gpu) && (use_gpu=="GPU" || use_gpu=="TRUE") )
	use_rgtsvm <- TRUE;
if (!is.na(use_gpu) && (use_gpu=="FALSE") )
	use_rgtsvm <- FALSE;

#src_path = as.character(args[7]);

cat("Bigwig(plus):", ps_plus_path, "\n");
cat("Bigwig(minus):", ps_minus_path, "\n");
cat("Output:", outfile, "\n");
cat("dREG model:", args[4], "\n");
cat("ncores:", ncores, "\n");
cat("GPU:", use_rgtsvm, "\n");
#cat("SRC PATH:", src_path, "\n");

## Now scan all positions in the genome ...
cat("1) -------- Checking the informative positions\n");
load(args[4]);

run.time <- system.time(r <- peak_calling( svm, gdm, ps_plus_path, ps_minus_path, ncores=ncores, use_rgtsvm=use_rgtsvm));

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

make_index_gz( r$infp_bed, out.file1 );
make_index_gz( r$peak_bed, out.file2 );
make_index_gz( r$peak_bed[,c(1:4)], out.file3 );
make_index_gz( r$peak_broad, out.file4 );


#system( paste("tar -cvzf ", outfile, ".dREG.tar.gz", " ", outfile, ".dREG.*", sep="") );
#cat("Result:", paste(outfile, ".dREG.tar.gz", sep=""), "\n");

cat("Done!\n");
