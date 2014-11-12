## Gets the mean phyloP scores inside of each functional class.
##

## Generate random:
## randomBed -l 500 -n 10000 -g hg19.chInfo | sort-bed - > hg19.random.bed

require(bigWig)
pp <- load.bigWig("/usr/data/hg19/phyloP_placMammals/merge/hg19.phyloP.placMammals.46way.bw")

beds <- c("gm12878.plus.dREG.bed", "gm12878.plus.dREG.PROONLY.bed", "gm12878.plus.dREG.ENHONLY.bed", 
	"gm12878.plus.DNAse.bed", "gm12878.chromHMM.only.bed", "gm12878.insulator.bed", "hg19.random.bed")

cols <- c("#e30000", "#8e0000", "#ff1c1c", "#000000", "#4f81bd", "#008e00", "gray")
step=10

## Define the mean phyloP score as the weighted mean by element size of the maxs.
maxs <- rep(0, length(beds))
vpl<- list()
mp <- list()
for(i in 1:length(beds)) {
 bed <- read.table(beds[i])
 bpq <- bed.region.bpQuery.bigWig(pp, bed, op= "max")
 vpl[[i]]<- bpq
 maxs[i] <- weighted.mean(x= bpq, w=(bed[,3]-bed[,2]))
 mp[[i]] <- meta.subsample(bed, pp, pp, 500, step, do.sum=TRUE)[[4]]/step
}
maxs

pdf("phyloP.placMammals.46way.pdf")
 plot(mp[[1]], type="l", col="#e30000", ylim=c(0.2,0.5))
 points(mp[[2]], type="l", col="#8e0000")
 points(mp[[3]], type="l", col="#ff1c1c")
 points(mp[[4]], type="l", col="#000000")
 points(mp[[5]], type="l", col="#4f81bd")
 points(mp[[6]], type="l", col="#008e00")
 points(mp[[7]], type="l", col="gray")

 library(vioplot)
 vioplot(vpl[[1]], vpl[[2]], vpl[[3]], vpl[[4]], vpl[[5]], vpl[[6]], vpl[[7]], col=cols)
dev.off()

