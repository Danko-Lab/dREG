setwd("/home/cgd24/storage/home/work/tss_detector/assayOverlap")

require(bigWig)

chrom <- read.table("k562.chromHMM.only.bed")[,1:3]
hs    <- read.table("k562.plus.DNAse.bed")[,1:3]
erna  <- read.table("k562.plus.dREG.bed")[,1:3]
erna_enh  <- read.table("k562.plus.dREG.ENHONLY.bed")[,1:3]
ins   <- read.table("k562.insulator.bed")[,1:3]  

doit <- function(file_name, stp=25, halfWindow= 5000, path= "/home/cgd24/storage/data/hg19/k562/histones/", ...) {
	bw <- load.bigWig(paste(path,file_name,sep=""))
	erna_meta <- metaprofile.bigWig(center.bed(erna, halfWindow, halfWindow), bw, step=stp)[[5]]/stp
	erna_enh_meta <- metaprofile.bigWig(center.bed(erna_enh, halfWindow, halfWindow), bw, step=stp)[[5]]/stp
	hs_meta <- metaprofile.bigWig(center.bed(hs, halfWindow, halfWindow), bw, step=stp)[[5]]/stp
	chrom_meta <- metaprofile.bigWig(center.bed(chrom, halfWindow, halfWindow), bw, step=stp)[[5]]/stp
	ins_meta <- metaprofile.bigWig(center.bed(ins, halfWindow, halfWindow), bw, step=stp)[[5]]/stp
	
    N = length(erna_meta)
    x = ((1:N) - N/2)* stp
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
 
 doit(file_name= "wgEncodeSydhNsomeK562Sig.bigWig", path="/home/cgd24/storage/data/hg19/k562/sydh_mnase/", main="MNase")
 doit(file_name= "wgEncodeSydhNsomeK562Sig.bigWig", halfWindow= 2000, path="/home/cgd24/storage/data/hg19/k562/sydh_mnase/", main="MNase")
dev.off()

