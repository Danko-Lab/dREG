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
chrom_info <- args[6]

ncores <- as.integer(args[7])
if (is.na(ncores)) ncores <- 1;

use_rgtsvm <- TRUE;
use_gpu <- toupper(as.character(args[8]))
if (!is.na(use_gpu) && use_gpu=="GPU")
	use_rgtsvm <- TRUE;

## Now scan all positions in the genome ...
cat("1) -------- Checking the informative positions\n");
inf_positions <- get_informative_positions(ps_plus_path, ps_minus_path, depth= 0, step=50, use_ANDOR=TRUE, use_OR=FALSE) ## Get informative positions.
cat("Genome Loci=", NROW(inf_positions), "\n");


cat("2) -------- Prediction in dREG model\n");
t <- system.time( pred_val<- eval_reg_svm(gdm, asvm, inf_positions, ps_plus_path, ps_minus_path, batch_size= 50000, ncores=ncores, use_rgtsvm=use_rgtsvm) )
cat("Running time [User]:", t[1], "[System]:", t[2], "[Elapsed]:", t[3], "\n");


cat("3) -------- write prediction data\n");
pred_data <- data.frame( inf_positions, pred_val )
file.dreg.pred.gz <- paste( outfile, ".dREG.pred.gz", sep="" );
gz1 <- gzfile( file.dreg.pred.gz, "w");
write.table(pred_data, gz1, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t");
close(gz1)

cat("4) -------- Merge peaks\n");
system( paste( "bash ", dirname(args[4]), "/writeBed.bsh", " 0.8 ", file.dreg.pred.gz, sep="") );

file.dreg.peak.gz <- paste(outfile, ".dREG.peak.gz", sep="");
file.rename(paste(file.dreg.pred.gz, ".bed.gz", sep=""), file.dreg.peak.gz)
unlink(file.dreg.pred.gz);

cat("5) -------- Refining in dREG-HD model\n");
require(dREG.HD);
## load drEG_HD model
load(args[5]);
t <- system.time( dREG_HD( bed_path = file.dreg.peak.gz, bigwig_plus = ps_plus_path, bigwig_minus = ps_minus_path, chromInfo=chrom_info, model=model, ncores = ncores, use_rgtsvm= use_rgtsvm ) )
cat("Running time [User]:", t[1], "[System]:", t[2], "[Elapsed]:", t[3], "\n");

file.rename( paste(file.dreg.peak.gz, "_imputedDnase.bw", sep=""), paste(outfile, "dREG.HD.imputedDnase.bw", sep=".") )
file.rename( paste(file.dreg.peak.gz, "_dREG_HD_relaxed.bed", sep=""), paste(outfile, "dREG.HD.relaxed.bed", sep=".") )
file.rename( paste(file.dreg.peak.gz, "_dREG_HD_stringent.bed", sep=""), paste(outfile, "dREG.HD.stringent.bed", sep=".") )

make_index_gz<-function(df_bed, outfile, file_id)
{
	write.table( df_bed, file=paste( outfile, ".", file_id, sep=""), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t");
	system(paste( "bgzip ", outfile, ".", file_id, sep="") );
	system(paste( "tabix -p bed ", outfile, ".", file_id, ".gz", sep="") );
}
make_index_gz(pred_data, outfile, "dREG.pred");

tb <- read.table(file.dreg.peak.gz, header = FALSE);
tb <- tb[,-4];

unlink( file.dreg.peak.gz );
make_index_gz( tb, outfile, "dREG.peak");

tb <- read.table(paste(outfile, "dREG.HD.relaxed.bed", sep="."));
tb <- data.frame(tb, 1);
make_index_gz( tb, outfile, "dREG.HD.relaxed.bed");

tb <- read.table(paste(outfile, "dREG.HD.stringent.bed", sep="."));
tb <- data.frame(tb, 1);
make_index_gz( tb, outfile, "dREG.HD.stringent.bed");

system( paste("tar -cvzf ", outfile, ".dREG.tar.gz", " ", outfile, ".dREG.*", sep="") );
cat("Result:", paste(outfile, ".dREG.tar.gz", sep=""), "\n");

