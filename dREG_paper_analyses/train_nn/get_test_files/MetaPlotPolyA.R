require(bigWig)

psplus  <- load.bigWig("/usr/data/GROseq.parser/hg19/k562/proseq/K562_unt.sort.bed.gz_plus.bw")
psminus <- load.bigWig("/usr/data/GROseq.parser/hg19/k562/proseq/K562_unt.sort.bed.gz_minus.bw")

k <- read.table("/usr/data/GROseq.parser/hg19/k562/polya_sites/k562.clusters.bed")
k <- k[k$V1 != "chrM",]

## Centered on 10kb upstream of the poly-A site.
l <- k
l$V3[k$V6 == "+"] <- l$V3[k$V6 == "+"]+10000
l$V2[k$V6 == "+"] <- l$V2[k$V6 == "+"]+10000
l$V3[k$V6 == "-"] <- l$V3[k$V6 == "-"]-10000
l$V2[k$V6 == "-"] <- l$V2[k$V6 == "-"]-10000

# Now choose the position from +5kb --> +15kb w/ highest difference in read densities in the upstream and downstream 5kb.
halfWindow= 10000
step= 500
kmany <- collect.many(l, psplus, psminus, halfWindow=halfWindow, step=step, do.sum=TRUE)

ncols <- 2*halfWindow/step
nwindows <- ncols/4

maxDiff <- sapply(1:NROW(k), function(x) {
	which.max(sapply(seq(ncols/4, 3*ncols/4, 1), function(y) {sum(kmany[x,(y-9):(y)])/(1+sum(kmany[x,(y+1):(y+10)]))}))
	})+10

## Get the position of the decrease.
m <- l
m$V3[k$V6 == "+"] <- k$V3[k$V6 == "+"]+step*maxDiff[k$V6 == "+"]
m$V2[k$V6 == "+"] <- k$V2[k$V6 == "+"]+step*maxDiff[k$V6 == "+"]
m$V3[k$V6 == "-"] <- k$V3[k$V6 == "-"]-step*maxDiff[k$V6 == "-"]
m$V2[k$V6 == "-"] <- k$V2[k$V6 == "-"]-step*maxDiff[k$V6 == "-"]

## Sanity check.
kcenters <- meta.subsample(m, psplus, psminus, halfWindow=halfWindow, step=step)
polyA    <- meta.subsample(k, psplus, psminus, halfWindow=halfWindow, step=step)

pdf("MetaPlots.polyA.pdf")
  par(mfrow=c(1,2))
  meta.plot(polyA, step=step, xlim=c(-10000,10000))
  meta.plot(kcenters, step=step, xlim=c(-10000,10000))
dev.off()

## Save yhose that dont overlap with tss

## Make a post-poly-A state.
n <- k
n$V3[k$V6 == "+"] <- m$V3[k$V6 == "+"]
n$V2[k$V6 == "-"] <- m$V2[k$V6 == "-"]

## get tss.
tss <- read.table("/home/cgd24/work/tss_detector/k562.predictions.bed.gz")

## intersect
require(rphast)
postpolya <- feat(seqname= n[,1], start= n[,2], end= n[,3])
tss <- feat(seqname= tss[,1], start= tss[,2], end= tss[,3])
nol <- overlap.feat(x= postpolya, filter= tss, overlapping = FALSE)
nool_indx <- match(paste(nol$seqname, nol$start, nol$end), paste(postpolya$seqname, postpolya$start, postpolya$end))

options("scipen"=100, "digits"=4)

## Write gene end state -- all post-polyA to max(dip) excluding regions w/ TSS in between.
write.table(m[nool_indx,], pipe(" sort-bed - > geneEnd.bed"), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

## Write 
write.table(n, pipe(" sort-bed - > postPolyA.bed"), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t") ## intersect this with the tss to get the post poly a state.


q("no")
################
## One last sanity check!
ge <- read.table("geneEnd.bed")
require(bigWig)

psplus  <- load.bigWig("/usr/data/GROseq.parser/hg19/k562/proseq/K562_unt.sort.bed.gz_plus.bw")
psminus <- load.bigWig("/usr/data/GROseq.parser/hg19/k562/proseq/K562_unt.sort.bed.gz_minus.bw")

halfWindow= 10000
step= 500

ge_c <- meta.subsample(ge, psplus, psminus, halfWindow=halfWindow, step=step)

meta.plot(ge_c, step=step, xlim=c(-10000,10000))

