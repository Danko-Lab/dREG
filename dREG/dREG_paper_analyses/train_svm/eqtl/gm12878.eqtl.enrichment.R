##
## Barplots of the number of eQTLs per detected site.
## 
## Compares between DNAse-1, Ernst' chromHMM, and ... 

#######################################################
## Get intersection between eQTL and each class.

eqtl <- read.table("eqtl.bed") ## Get eQTLs.
tss <- read.table("gm12878.bed")##gm12878.tss.10fdr.bed") ## Get TSS ... ##### 
dnase1 <- read.table("dnase.narrowpeak.bed") ## Get DNAse-1. 
chromhmm <- read.table("chromHMM.bed") ## Get chromHMM.
tmp.bed <- read.table("tmp.bed")

getOverlap <- function(dataset, eqtl) {
  require(rphast)
  eqtl_feat  <- feat(seqname= eqtl[,1], start= eqtl[,2], end= eqtl[,3])
  data_feat <- feat(seqname= dataset[,1], start= dataset[,2], end= dataset[,3])
  ol <- overlap.feat(x= eqtl_feat, filter= data_feat)
  pos_indx <- match(paste(ol$seqname, ol$start, ol$end), paste(eqtl_feat$seqname, eqtl_feat$start, eqtl_feat$end))
  return(NROW(pos_indx))
}

sitesRecovered <- list("GROseq TSS"= getOverlap(tss, eqtl), "Dnase-1"= getOverlap(dnase1, eqtl), "ChromHMM"= getOverlap(chromhmm, eqtl))#, "chromHMM.ActiveProm.Or.StrongEnh"= getOverlap(tmp.bed, eqtl))
n_sites <- list(tss= NROW(tss), dnase1= NROW(dnase1), chromhmm= NROW(chromhmm))#, chromHMM2= NROW(tmp.bed))

data_df <- data.frame(Assay= names(sitesRecovered), SitesRecovered= as.double(sitesRecovered)/as.double(n_sites)) ## Plot ...
data_df$Assay <- factor(data_df$Assay, levels= levels(data_df$Assay)[c(3,2,1)])
data_df

# require(ggplot2)
# library(reshape2)
# http://learnr.wordpress.com/2009/03/17/ggplot2-barplots/
# http://stackoverflow.com/questions/10352894/barplot-using-ggplot2
#a <- ggplot(data_df, aes(Type, Changed)) + xlab("") + ylab("Fraction of Transcripts Changed") +
#	geom_bar(aes(fill=Type), stat="identity") + 
#	scale_fill_brewer(palette = "Set1")

pdf("fig5a.pdf")

require(lattice)
barchart(SitesRecovered~Assay, data_df, auto.key = list(columns = 1), 
	xlab= list(label= "Method", cex=2.05), 
	ylab= list(label="eQTLs per site", cex=2.05), 
	lwd=1, pch=1, scales=list(x=list(cex=1.95), y=list(cex=1.95)),
	col=c("#A25540","#95BD55","#7F9D9B","#9D5DAB"))	

dev.off()

#######################################################
## Plot recovery as a function of threshold.
# th <- seq(60, 95, 5)
# recovery <- list()
# sites <- list()
# for(i in th) {
  # tssi <- read.table(paste("gm12878.tss",i,".bed",sep=""))
  # recovery[[paste(i)]] <- getOverlap(tssi, eqtl)
  # sites[[paste(i)]] <- NROW(tssi)
# }
# tss_eqtl_frac <- as.double(recovery)/as.double(sites)

# pdf("figS5ab.pdf")

# par(mfrow=c(1,2), cex=1.2, cex.lab=1.2, cex.axis=1.05, lwd=2, font.axis=2, font.lab=2, font=2)
# plot(th, tss_eqtl_frac, type="b", main="eQTL per Site", xlab="TSS Threshold", ylab="eQTLs per site",
	# ylim=c(min(data_df$SitesRecovered, tss_eqtl_frac),max(data_df$SitesRecovered, tss_eqtl_frac)), 
	# col="#A25540", pch=19)
# abline(h=(data_df$SitesRecovered[data_df$Assay == "Dnase-1"]), lty="dashed", col="#95BD55")
# abline(h=data_df$SitesRecovered[data_df$Assay == "ChromHMM"], lty="dashed", col="#7F9D9B")

# plot(th, as.double(recovery), type="b", main="eQTL total recovery", xlab="TSS Threshold", ylab="Total eQTLs recovered",
	# ylim=c(min(as.double(sitesRecovered), as.double(recovery)),max(as.double(sitesRecovered), as.double(recovery))), 
	# col="#A25540", pch=19)
# abline(h=(sitesRecovered[["Dnase-1"]]), lty="dashed", col="#95BD55")
# abline(h=(sitesRecovered[["ChromHMM"]]), lty="dashed", col="#7F9D9B")

# dev.off()

#######################################################
## Plot the density of eQTLs near sites...

require(bigWig)
stepsize <- 200
halfWindow <- 5000
x_axis <- seq(-halfWindow+stepsize/2,halfWindow-stepsize/2,stepsize)

eqbw <- load.bigWig("eqtl.bw")
eq_meta_tss<- meta.accum(bed= data.frame(tss[,c(1:3)], "+"), bigWig.plus= eqbw, halfWindow=halfWindow, step=stepsize, do.sum = TRUE)
eq_meta_dnase<- meta.accum(bed= dnase1[,c(1:3,6)], bigWig.plus= eqbw, halfWindow=halfWindow, step=stepsize, do.sum = TRUE)
eq_meta_chromhmm<- meta.accum(bed= chromhmm[,c(1:3,6)], bigWig.plus= eqbw, halfWindow=halfWindow, step=stepsize, do.sum = TRUE)
eq_meta_pr<- meta.accum(bed= tmp.bed[,c(1:3,6)], bigWig.plus= eqbw, halfWindow=halfWindow, step=stepsize, do.sum = TRUE)
ylim=c(min(eq_meta_tss,eq_meta_dnase,eq_meta_chromhmm,eq_meta_pr), max(eq_meta_tss,eq_meta_dnase,eq_meta_chromhmm,eq_meta_pr))

pdf("fig5b.pdf")
 par(cex=1.2, cex.lab=1.2, cex.axis=1.15, lwd=3)
 plot(seq(-halfWindow+stepsize/2,halfWindow-stepsize/2,stepsize), eq_meta_tss, type="l", main="eQTLs on GRO-seq TSS", col="#A25540", ylab="eQTL density", xlab="Distance from center [bp]", ylim=ylim)
 points(seq(-halfWindow+stepsize/2,halfWindow-stepsize/2,stepsize), eq_meta_dnase, type="l", col="#95BD55")
 points(seq(-halfWindow+stepsize/2,halfWindow-stepsize/2,stepsize), eq_meta_chromhmm, type="l", col="#7F9D9B")
#points(seq(-halfWindow+stepsize/2,halfWindow-stepsize/2,stepsize), eq_meta_pr, type="l", col="#7F9D9B")
 legend(x= 1500, y= 0.0025, legend= c("GROseq", "DNAse-1", "chromHMM"), col= c("#A25540", "#95BD55", "#7F9D9B"), lwd= c(2, 2, 2))
dev.off()

## And for Andre, do from the center of GRO-cap double peaks.
get_grocap_pair <- function(plus, minus) {
 grocap_pair_p <- read.table(plus)
 grocap_pair_m <- read.table(minus)
 grocap_pair <- grocap_pair_p
 grocap_pair[,2] <- grocap_pair_m[,2]
 return(grocap_pair)
}

grocap_pair <- get_grocap_pair(plus= "andre_hmm/hg19.gm12878.new_hmm2b.post2.pair_plus.bed", minus= "andre_hmm/hg19.gm12878.new_hmm2b.post2.pair_minus.bed")
eq_meta_grocap<- meta.accum(bed= grocap_pair[,c(1:3,6)], bigWig.plus= eqbw, halfWindow=halfWindow, step=stepsize)

pdf("fig5b_wGROcap.pdf")
par(cex=1.2, cex.lab=1.2, cex.axis=1.15, lwd=3)

plot(x_axis, eq_meta_tss, type="l", main="eQTL density near regulatory sites", col="#A25540", ylab="eQTL density", xlab="Distance from center [bp]", 
	ylim=c(min(eq_meta_grocap), max(eq_meta_grocap)))
points(x_axis, eq_meta_dnase, type="l", col="#95BD55")
points(x_axis, eq_meta_chromhmm, type="l", col="#7F9D9B")
points(x_axis, eq_meta_grocap, type="l", col="red")
legend(x= 1500, y= 0.004, legend= c("GROseq", "GROcap", "DNAse-1", "chromHMM"), col= c("#A25540", "red", "#95BD55", "#7F9D9B"), lwd= c(2, 2, 2, 2))
dev.off()

## And look near stability classes.
stepsize <- 100
halfWindow <- 1000
x_axis <- seq(-halfWindow+stepsize/2,halfWindow-stepsize/2,stepsize)

SS <- get_grocap_pair(plus= "andre_hmm/hg19.gm12878.new_hmm2b.post2.SS_plus.bed", minus= "andre_hmm/hg19.gm12878.new_hmm2b.post2.SS_minus.bed")
SU <- get_grocap_pair(plus= "andre_hmm/hg19.gm12878.new_hmm2b.post2.SU_plus.bed", minus= "andre_hmm/hg19.gm12878.new_hmm2b.post2.SU_minus.bed")
US <- get_grocap_pair(plus= "andre_hmm/hg19.gm12878.new_hmm2b.post2.US_plus.bed", minus="andre_hmm/hg19.gm12878.new_hmm2b.post2.US_minus.bed")
UU <- get_grocap_pair(plus= "andre_hmm/hg19.gm12878.new_hmm2b.post2.UU_plus.bed", minus="andre_hmm/hg19.gm12878.new_hmm2b.post2.UU_minus.bed")

meta_SS<- meta.accum(bed= SS[,c(1:3,6)], bigWig.plus= eqbw, halfWindow=halfWindow, step=stepsize, do.sum = TRUE)
meta_SU<- meta.accum(bed= SU[,c(1:3,6)], bigWig.plus= eqbw, halfWindow=halfWindow, step=stepsize, do.sum = TRUE)
meta_US<- meta.accum(bed= US[,c(1:3,6)], bigWig.plus= eqbw, halfWindow=halfWindow, step=stepsize, do.sum = TRUE)
meta_UU<- meta.accum(bed= UU[,c(1:3,6)], bigWig.plus= eqbw, halfWindow=halfWindow, step=stepsize, do.sum = TRUE)

pdf("eqtl.grocap_SS.US.UU.pdf")
par(cex=1.2, cex.lab=1.2, cex.axis=1.15, lwd=3)

plot(x_axis, meta_SS, type="l", main="eQTLs on GRO-cap pairs", col="red", ylab="eQTL density", xlab="Distance from center [bp]", 
	ylim=c(min(meta_SS, meta_SU, meta_US, meta_UU), max(meta_SS, meta_SU, meta_US, meta_UU)))
points(x_axis, (rev(meta_SU)+meta_US)/2, type="l", col="green")
points(x_axis, meta_UU, type="l", col="blue")
abline(v=0, lty="dotted", col="dark gray")
legend(x= 500, y= 0.005, legend= c("SS", "US", "UU"), col= c("red", "green", "blue"), lwd= c(2, 2, 2))

dev.off()


#######################################################
## Finally break out eqtl by the 'overlaps' defined in fig. 3...

setwd("/home/cgd24/work/tss_detector/assayOverlap")

prefix <- "gm12878"
chrom <- read.table(paste(prefix, ".chromHMM.only.bed.gz", sep=""))
dnase <- read.table(paste(prefix, ".plus.DNAse.bed.gz", sep=""))
dreg  <- read.table(paste(prefix, ".plus.dREG.bed.gz", sep=""))
dregE <- read.table(paste(prefix, ".plus.dREG.ENHONLY.bed.gz", sep=""))

require(bigWig)
stepsize <- 200
halfWindow <- 5000
x_axis <- seq(-halfWindow+stepsize/2,halfWindow-stepsize/2,stepsize)

eqbw <- load.bigWig("eqtl.bw")
m_chrom<- meta.accum(bed= chrom, bigWig.plus= eqbw, halfWindow=halfWindow, step=stepsize, do.sum = TRUE)
m_dnase<- meta.accum(bed= dnase, bigWig.plus= eqbw, halfWindow=halfWindow, step=stepsize, do.sum = TRUE)
m_dreg <- meta.accum(bed= dreg, bigWig.plus= eqbw, halfWindow=halfWindow, step=stepsize, do.sum = TRUE)
m_dregE<- meta.accum(bed= dregE, bigWig.plus= eqbw, halfWindow=halfWindow, step=stepsize, do.sum = TRUE)

ylim=c(min(m_chrom,m_dnase,m_dreg), max(m_chrom,m_dnase,m_dreg))
xlim=c(-5000,5000)

pdf("fig4c.eqtlEnrichmentByClass.pdf")
 par(cex=1.2, cex.lab=1.2, cex.axis=1.15, lwd=3)
 plot(seq(-halfWindow+stepsize/2,halfWindow-stepsize/2,stepsize), m_dreg, type="l", main="eQTLs on GRO-seq TSS", col="#e30000", ylab="eQTL density", xlab="Distance from center [bp]", ylim=ylim, xlim=xlim)
 points(seq(-halfWindow+stepsize/2,halfWindow-stepsize/2,stepsize), m_dregE, type="l", col="red")
 points(seq(-halfWindow+stepsize/2,halfWindow-stepsize/2,stepsize), m_dnase, type="l", col="#000000")
 points(seq(-halfWindow+stepsize/2,halfWindow-stepsize/2,stepsize), m_chrom, type="l", col="#4f81bd")
 legend(x= 1500, y= 0.0025, legend= rev(c("+dREG", "+DNAse", "chromHMM")), col= rev(c("#e30000", "#000000", "#4f81bd")), lwd= c(2, 2, 2))
dev.off()
