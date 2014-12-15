library(lattice)
setwd("/home/cgd24/storage/home/work/tss_detector/assayOverlap")

## Fig. 3C.  Compare chromHMM, DNAse1, gm12878 with regards to TF binding.

prefix <- "k562"
chrom <- read.table(paste(prefix, ".chromhmm.tsv.gz", sep=""))
dnase <- read.table(paste(prefix, ".dnase.tsv.gz", sep=""))
ins <- read.table(paste(prefix, ".insulator.tsv.gz", sep=""))
dreg  <- read.table(paste(prefix, ".dreg.tsv.gz", sep=""))
dregE <- read.table(paste(prefix, ".dreg.ENHONLY.tsv.gz", sep=""))

## Remove columns corresponding to GTFs, Pol II, TBP
col_names <- read.table(paste(prefix,".colNames.txt", sep=""))$V3

## Collapse same TFs by taking 1 if any binding occurs in all peak calls.
chrom_combined <- NULL
dnase_combined <- NULL
ins_combined <- NULL
dreg_combined <- NULL
dregE_combined <- NULL
for(i in unique(col_names[grep("POL|GTF|TBP", col_names, invert=TRUE)])) {
  if(sum(col_names == i)>1) {
    chrom_combined <- cbind(chrom_combined, as.integer(rowSums(chrom[,3+which(col_names == i)])>0))
    dnase_combined <- cbind(dnase_combined, as.integer(rowSums(dnase[,6+which(col_names == i)])>0))
    ins_combined <- cbind(ins_combined, as.integer(rowSums(ins[,3+which(col_names == i)])>0))
    dreg_combined <- cbind(dreg_combined, as.integer(rowSums(dreg[,5+which(col_names == i)])>0))
    dregE_combined <- cbind(dregE_combined, as.integer(rowSums(dregE[,5+which(col_names == i)])>0))
  }
  else {
    chrom_combined <- cbind(chrom_combined, chrom[,3+which(col_names == i)])
    dnase_combined <- cbind(dnase_combined, dnase[,6+which(col_names == i)])
    ins_combined <- cbind(ins_combined, ins[,3+which(col_names == i)])
    dreg_combined <- cbind(dreg_combined, dreg[,5+which(col_names == i)])
    dregE_combined <- cbind(dregE_combined, dregE[,5+which(col_names == i)])
  }
}

chrom <- rowSums(chrom_combined)
dnase <- rowSums(dnase_combined)
ins   <- rowSums(ins_combined)
dreg  <- rowSums(dreg_combined)
dregE <- rowSums(dregE_combined)

pdf("Fig3C.nTFs.pdf")

par(mfrow=c(5,1))
hist(chrom, 50, main="chrom")
hist(dnase, 50, main="dnase")
hist(ins, 50, main="ins")
hist(dregE, 50, main="dregE")
hist(dreg, 50, main="dreg")

df <- data.frame(groups= c("0", "1", "2", "3-5", "6-10", "11-20", "21-30", ">30"), 
                 chrom= c(sum(chrom==0), sum(chrom==1), sum(chrom==2), sum(chrom>=3 & chrom<=5), sum(chrom>=6 & chrom<=10), sum(chrom>=11 & chrom<=20), sum(chrom>=21 & chrom<=30), sum(chrom>30))/NROW(chrom),
				 dnase= c(sum(dnase==0), sum(dnase==1), sum(dnase==2), sum(dnase>=3 & dnase<=5), sum(dnase>=6 & dnase<=10), sum(dnase>=11 & dnase<=20), sum(dnase>=21 & dnase<=30), sum(dnase>30))/NROW(dnase),
				 ins= c(sum(ins==0), sum(ins==1), sum(ins==2), sum(ins>=3 & ins<=5), sum(ins>=6 & ins<=10), sum(ins>=11 & ins<=20), sum(ins>=21 & ins<=30), sum(ins>30))/NROW(ins),
				 dregE=  c(sum(dregE==0), sum(dregE==1), sum(dregE==2), sum(dregE>=3 & dregE<=5), sum(dregE>=6 & dregE<=10), sum(dregE>=11 & dregE<=20), sum(dregE>=21 & dregE<=30), sum(dregE>30))/NROW(dregE),
				 dreg=  c(sum(dreg==0), sum(dreg==1), sum(dreg==2), sum(dreg>=3 & dreg<=5), sum(dreg>=6 & dreg<=10), sum(dreg>=11 & dreg<=20), sum(dreg>=21 & dreg<=30), sum(dreg>30))/NROW(dreg)
				 )
df

daf2 <- data.frame(ntfs= rep(as.factor(c(0:26)), 3), type= as.factor(c(rep("chrom", 27), rep("dnase", 27), rep("dreg", 27))),
			nTFs= c(c(sapply(c(0:25), function(x) {sum(chrom == x)}), sum(chrom>25))/ NROW(chrom), 
			c(sapply(c(0:25), function(x) {sum(dnase == x)}), sum(dnase>25))/ NROW(dnase), 
			c(sapply(c(0:25), function(x) {sum(dregE == x)}), sum(dregE>25))/ NROW(dregE)))
daf2

daf <- data.frame(ntfs= rep(as.factor(c(0:21)), 4), type= as.factor(c(rep("chrom", 22), rep("dnase", 22), rep("ins", 22), rep("dregE", 22))),
			nTFs= c(c(sapply(c(0:20), function(x) {sum(chrom == x)}), sum(chrom>20))/ NROW(chrom), 
			c(sapply(c(0:20), function(x) {sum(dnase == x)}), sum(dnase>20))/ NROW(dnase), 
			c(sapply(c(0:20), function(x) {sum(ins == x)}), sum(ins>20))/ NROW(ins), 
			c(sapply(c(0:20), function(x) {sum(dregE == x)}), sum(dregE>20))/ NROW(dregE)))
daf

dafp <- data.frame(ntfs= rep(as.factor(c(0:21)), 4), type= as.factor(c(rep("chrom", 22), rep("dnase", 22), rep("ins", 22), rep("dreg", 22))),
			nTFs= c(c(sapply(c(0:20), function(x) {sum(chrom == x)}), sum(chrom>20))/ NROW(chrom), 
			c(sapply(c(0:20), function(x) {sum(dnase == x)}), sum(dnase>20))/ NROW(dnase), 
			c(sapply(c(0:20), function(x) {sum(ins == x)}), sum(ins>20))/ NROW(ins), 
			c(sapply(c(0:20), function(x) {sum(dreg == x)}), sum(dreg>20))/ NROW(dreg)))
dafp


library(RColorBrewer)


my.settings <- list(
  superpose.polygon=list(col=c("#4f81bd", "#000000", "#e30000", "#218a15"), border="transparent"),
  strip.background=list(col=c("#4f81bd", "#000000", "#e30000", "#218a15")),
  strip.border=list(col="black")
)


barchart(nTFs~ntfs, groups=type, data=daf, ylim=c(0,0.7), ylab="Fraction of Peaks", xlab="Number of TFs", 
		par.settings = my.settings, par.strip.text=list(col="white", font=2),
		auto.key=list(space="top", columns=4))

barchart(nTFs~ntfs, groups=type, data=dafp, ylim=c(0,0.7), ylab="Fraction of Peaks", xlab="Number of TFs", 
		par.settings = my.settings, par.strip.text=list(col="white", font=2),
		auto.key=list(space="top", columns=4))

barchart(nTFs~ntfs, groups=type, data=daf2, ylim=c(0,0.7), ylab="Fraction of Peaks", xlab="Number of TFs", 
		par.settings = my.settings, par.strip.text=list(col="white", font=2),
		auto.key=list(space="top", columns=3))

dev.off()

