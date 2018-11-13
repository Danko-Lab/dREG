library(bigWig);

library("tools")
options("scipen"=100, "digits"=4)
options(width=250)

bed_intersect <- function( score_bed, refer_bed ){
	file.refer <- tempfile(fileext=".bed");
	file.score <- tempfile(fileext=".bed");
	write.table( refer_bed[,c(1:3)], file=file.refer, quote=F, col.names=F, row.names=F, sep="\t");
	write.table( score_bed, file=file.score, quote=F, col.names=F, row.names=F, sep="\t");
	ret <- unique( read.table( pipe(paste("bedtools intersect -a ", file.score,  " -b ", file.refer, " -wa ", sep=" ")) ) );
	unlink(file.refer);
	unlink(file.score);
	return(ret);
}

bed_subtract <- function( score_bed, refer_bed ){
	file.refer <- tempfile(fileext=".bed");
	file.score <- tempfile(fileext=".bed");
	write.table( refer_bed[,c(1:3)], file=file.refer, quote=F, col.names=F, row.names=F, sep="\t");
	write.table( score_bed, file=file.score, quote=F, col.names=F, row.names=F, sep="\t");
	ret <- unique( read.table( pipe(paste("bedtools intersect -a ", file.score,  " -b ", file.refer, " -v ", sep=" ")) ) );
	unlink(file.refer);
	unlink(file.score);
	return(ret);
}


file.bw.G1 <- c( "../k562/K562_unt.sort.bed.gz_plus.bw", "../k562/K562_unt.sort.bed.gz_minus.bw");
file.G1.plus <- file.bw.G1[1];
file.G1.minus <- file.bw.G1[2];

file.K562.DHS.peak  <- "./k562.merged.broad.peak.bed";
file.grocap="../k562/hg19.k562.new_hmm2b.post2.bed"
TRE.bed <- read.table("/local/workdir/zw355/proj/prj10-dreg/new-rf-201803/G1/G1.dREG.peak.score.bed.gz");

DHS.bed <- read.table(file.K562.DHS.peak)
NONDHS.bed <- bed_subtract(TRE.bed, DHS.bed);
dREG.DHSplus.bed <- bed_intersect(TRE.bed, DHS.bed);
dREG.DHSminus.bed <- NONDHS.bed;

dREG.DHSplus.bed[,3] <- dREG.DHSplus.bed[,3] + 250
dREG.DHSplus.bed[,2] <- dREG.DHSplus.bed[,2] - 250
dREG.DHSminus.bed[,3] <- dREG.DHSminus.bed[,3] + 250
dREG.DHSminus.bed[,2] <- dREG.DHSminus.bed[,2] - 250

bw.plus <- load.bigWig( file.G1.plus);
bw.minus <- load.bigWig( file.G1.minus);

reads.plus <- bed.region.bpQuery.bigWig( bw.plus, dREG.DHSplus.bed[,c(1:3)], op = "sum", abs.value = TRUE)
             +bed.region.bpQuery.bigWig( bw.minus,dREG.DHSplus.bed[,c(1:3)], op = "sum", abs.value = TRUE);    
reads.minus<- bed.region.bpQuery.bigWig( bw.plus, dREG.DHSminus.bed[,c(1:3)], op = "sum", abs.value = TRUE)
             +bed.region.bpQuery.bigWig( bw.minus,dREG.DHSminus.bed[,c(1:3)], op = "sum", abs.value = TRUE);

unload.bigWig(bw.plus);
unload.bigWig(bw.minus);

reads <- data.frame(count=c(reads.plus, reads.minus)+1, type=c(rep("dREG+DHS+", NROW(reads.plus)), rep("dREG+DHS-", NROW(reads.minus)) ) );

pdf("fig-S7.pdf")
boxplot(count~type,data=reads, main="",   xlab="", ylab="Read counts", log="y" ) 
dev.off();

