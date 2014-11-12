## Get ratios of H3K4me3 to H3K4me1.

args <- commandArgs(trailingOnly=TRUE)
bed_p <- args[1]
me3_p <- args[2]
me1_p <- args[3]

require(bigWig)
bed <- read.table(bed_p)[,c(1:3)]

me3 <- load.bigWig(me3_p)
me1 <- load.bigWig(me1_p)

c3 <- bedQuery.bigWig(bed, me3, gapValue=0)
c1 <- bedQuery.bigWig(bed, me1, gapValue=0)

countReads <- function(bw) {
  t_sum=0
  for(i in 1:NROW(bw$chroms)) {
    t_sum<- t_sum+chromStepSum.bigWig(bw, bw$chroms[i], bw$chromSizes[i], 0)
  }
  return(t_sum)
}

nf3<-1# countReads(me3)
nf1<-1# countReads(me1)

pc <- 0.001
lr <- log(((c3+pc)/nf3)/((c1+pc)/nf1))
summary(lr)

#hist(lr, 150)
write.table(lr, "lr.tmp", quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t")

