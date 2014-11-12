

##
ss <- read.table("superset.nums.bed")

## Using these values to construct a Venn diagriam by hand.

## dREG is ss[,4]
sum(ss[,4] == 1 & ss[,5] == 1 & ss[,6] == 1) # 22254
sum(ss[,4] == 1 & ss[,5] == 0 & ss[,6] == 0) # 6784
sum(ss[,4] == 1 & ss[,5] == 1 & ss[,6] == 0) # 645
sum(ss[,4] == 1 & ss[,5] == 0 & ss[,6] == 1) # 3068

## DNASE-1 is ss[,5]
sum(ss[,4] == 1 & ss[,5] == 1 & ss[,6] == 1) # 22254
sum(ss[,4] == 0 & ss[,5] == 1 & ss[,6] == 0) # 47963
#sum(ss[,4] == 1 & ss[,5] == 1 & ss[,6] == 0) # 645 ## REPEAT.
sum(ss[,4] == 0 & ss[,5] == 1 & ss[,6] == 1) # 23287

## chromHMM is ss[,6]
#sum(ss[,4] == 1 & ss[,5] == 1 & ss[,6] == 1) # 22254 ## REPEAT.
sum(ss[,4] == 0 & ss[,5] == 0 & ss[,6] == 1) # 65881
#sum(ss[,4] == 0 & ss[,5] == 1 & ss[,6] == 1) # 23287 ## REPEAT.
#sum(ss[,4] == 1 & ss[,5] == 0 & ss[,6] == 1) # 3068 ## REPEAT.

## H3k27ac is ss[,7]
##
## This is interesting ... does this mean that transcription preceeds H3K27ac?
## Hmm.  h3k27ac covers more of hg19 (2.685% [h3k12ac] vs. 1.283% [dreg] vs. 1.041% [intersection]).  

#sh <- read.table("superset.nums.hist.bed")
#sum(sh[,4] == 1 & sh[,5] == 1) # 18817 (48% of dREG, 72% of H3k27ac).
#sum(sh[,4] == 1 & sh[,5] == 0) # 20178
#sum(sh[,4] == 0 & sh[,5] == 1) # 7329

######################################################
## Print overlaps for Fig. 3B.
library(lattice)

x_cex=y_cex=1.5
x_lab_cex=y_lab_cex=1.55
col=c("dark green", "dark blue")

overlap <- read.csv(text="Mark, Overlap
DNAse-1, 27
chromHMM, 22
H3K27ac, 72
H3K9ac, 86", header=TRUE)

barchart(Overlap~Mark,data=overlap, col="dark blue",
         scales=list(x=list(rot=65,cex=x_cex), y=list(rot=0, cex=y_cex)), ylim=c(0,100), ylab="Fraction of Sites Transcribed",
		 auto.key=list(columns = 2, space = "top"),
		 lwd=1, pch=1, xlab= list(cex=x_lab_cex))
