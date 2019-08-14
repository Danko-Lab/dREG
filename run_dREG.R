require(dREG)
options("scipen"=100, "digits"=4)

make_index_gz<-function( df_bed, out_file)
{
	file.tmp <- tempfile(fileext=".bed");
	write.table( df_bed, file=file.tmp, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t");
	system(paste( "sort-bed ", file.tmp, " | bgzip -f > ", out_file, ".gz", sep="") );
	system(paste( "tabix -f -p bed ", out_file, ".gz", sep="") );
	unlink(file.tmp)
}

make_bw<-function( df_bed, out_file)
{
    ## decrease the size
    df_bed[,4] <- round(df_bed[,4], 2);

    file.tmp <- tempfile(fileext=".bed");
    write.table( df_bed, file=file.tmp, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t");

    file.tmp2 <- tempfile(fileext=".bed");
    system(paste( "sort-bed ", file.tmp, ">", file.tmp2) );
    system(paste( "bedGraphToBigWig ", file.tmp2, out.chrom.info, out_file,  sep=" ") );

    unlink(file.tmp)
    unlink(file.tmp2)
}

make_chrom_info <- function( file.bw.plus, file.bw.minus, file.chrom.info)
{
    bw.plus <- load.bigWig(file.bw.plus);
    bw.minus <- load.bigWig(file.bw.minus);

    chrom <- rbind( cbind( bw.plus$chroms, bw.plus$chromSizes), cbind( bw.minus$chroms, bw.minus$chromSizes) );
    chr.size <- unlist( lapply( unique(chrom[,1]), function(chr){max( as.numeric( chrom[which(chrom[,1]==chr),2])) } ) );

    df.bed <- data.frame( V1=unique(chrom[,1]), V2=chr.size );
    write.table( df.bed, file=file.chrom.info, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t");
}



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

out.rdata<- paste(outfile, "dREG.rdata", sep=".")
## RDATA file is used for debug
## save(r, file=out.rdata );

out.raw.peak.bed<- paste(outfile, "dREG.raw.peak.bed", sep=".")
out.infp.bed <- paste(outfile, "dREG.infp.bed", sep=".")
out.infp.bw  <- paste(outfile, "dREG.infp.bw", sep=".")
out.score.bed<- paste(outfile, "dREG.peak.score.bed", sep=".")
out.score.bw <- paste(outfile, "dREG.peak.score.bw", sep=".")
out.prob.bed <- paste(outfile, "dREG.peak.prob.bed", sep=".")
out.prob.bw  <- paste(outfile, "dREG.peak.prob.bw", sep=".")
out.full.bed <- paste(outfile, "dREG.peak.full.bed", sep=".")
out.chrom.info <- paste(outfile, "chrom.info", sep=".")

r$infp_bed[,4] <- round( r$infp_bed[,4], 3 );
r$peak_bed[,4] <- round( r$peak_bed[,4], 3 );

make_chrom_info(ps_plus_path, ps_minus_path, out.chrom.info);
make_index_gz( r$raw_peak[ r$raw_peak$prob != -1, ], out.raw.peak.bed);
make_index_gz( r$infp_bed, out.infp.bed );
make_index_gz( r$peak_bed, out.full.bed );
make_index_gz( r$peak_bed[,c(1:4)], out.score.bed );

r$peak_bed[,5] <- 1 - r$peak_bed[,5];
make_index_gz( r$peak_bed[,c(1,2,3,5)], out.prob.bed );

make_bw( r$infp_bed[,1:4], out.infp.bw );
make_bw( r$peak_bed[,c(1:4)], out.score.bw );
make_bw( r$peak_bed[,c(1,2,3,5)], out.prob.bw );

#system( paste("tar -cvzf ", outfile, ".dREG.tar.gz", " ", outfile, ".dREG.*", sep="") );
#cat("Result:", paste(outfile, ".dREG.tar.gz", sep=""), "\n");

cat("Done!\n");