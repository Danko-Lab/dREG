setwd("/home/cgd24/work/tss_detector/assayOverlap")

chrom <- read.table("k562.chromHMM.only.bed")
hs    <- read.table("k562.plus.DNAse.bed")
erna  <- read.table("k562.plus.dREG.bed")
erna_enh  <- read.table("k562.plus.dREG.ENHONLY.bed")
ins   <- read.table("k562.insulator.bed")

require(bigWig)

doit <- function(file_name, halfWindow=5000, step=25, path= "/usr/data/GROseq.parser/hg19/k562/histones/", ...) {
	bw <- load.bigWig(paste(path,file_name,sep=""))
	erna_meta <- meta.subsample(erna, bw, bw, halfWindow=halfWindow, step=step, do.sum=TRUE)[[4]]/step
	erna_enh_meta <- meta.subsample(erna_enh, bw, bw, halfWindow=halfWindow, step=step, do.sum=TRUE)[[4]]/step
	hs_meta <- meta.subsample(hs, bw, bw, halfWindow=halfWindow, step=step, do.sum=TRUE)[[4]]/step
	chrom_meta <- meta.subsample(chrom, bw, bw, halfWindow=halfWindow, step=step, do.sum=TRUE)[[4]]/step
	ins_meta <- meta.subsample(ins, bw, bw, halfWindow=halfWindow, step=step, do.sum=TRUE)[[4]]/step
	
    N = length(erna_meta)
    x = ((1:N) - N/2)* step
	ylim=c(min(c(erna_meta, hs_meta, chrom_meta, ins_meta, erna_enh_meta)), max(c(erna_meta, hs_meta, chrom_meta, ins_meta, erna_enh_meta)))
	
	plot(x, erna_meta, type="l", ylim=ylim, col="#e30000", ...)
	points(x, hs_meta, type="l", ylim=ylim, col="#000000")
	points(x, chrom_meta, type="l", ylim=ylim, col="#4f81bd")
	points(x, ins_meta, type="l", ylim=ylim, col="#218a15")
	points(x, erna_enh_meta, type="l", ylim=ylim, col="grey")
}

pdf("histoneMetaPlots.overlap.pdf")
 doit(file_name= "wgEncodeBroadHistoneK562P300StdSig.bigWig", main="P300")
 doit(file_name= "wgEncodeBroadHistoneK562H3k27acStdSig.bigWig", main="H3k27ac")
 doit(file_name= "wgEncodeBroadHistoneK562H3k9acStdSig.bigWig", main="H3k9ac")

 doit(file_name= "wgEncodeBroadHistoneK562H3k4me1StdSig.bigWig", main="H3k4me1")
 doit(file_name= "wgEncodeBroadHistoneK562H3k4me3StdSig.bigWig", main="H3k4me3")

 doit(file_name= "wgEncodeBroadHistoneK562H3k27me3StdSig.bigWig", main="H3k27me3")
 
 doit(file_name= "wgEncodeBroadHistoneK562H2azStdSig.bigWig", main="H2az")
 
 doit(file_name= "wgEncodeSydhNsomeK562Sig.bigWig", path="/usr/data/GROseq.parser/hg19/k562/sydh_mnase/", main="MNase")
 doit(file_name= "wgEncodeSydhNsomeK562Sig.bigWig", halfWindow= 2000, path="/usr/data/GROseq.parser/hg19/k562/sydh_mnase/", main="MNase")
dev.off()

