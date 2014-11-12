require(bigWig)

psplus  <- load.bigWig("/usr/data/GROseq.parser/hg19/k562/proseq/K562_unt.sort.bed.gz_plus.bw")
psminus <- load.bigWig("/usr/data/GROseq.parser/hg19/k562/proseq/K562_unt.sort.bed.gz_minus.bw")

txn <- read.table("ernst_txn.bed")

## Get the expression of GRO-cap in each cluster.  
plus <- bedQuery.bigWig(txn, psplus, gapValue=0)*1000/(txn[,3]-txn[,2])
minus <- bedQuery.bigWig(txn, psminus, gapValue=0)*1000/(txn[,3]-txn[,2])

 rowMins <- function(x, y) {sapply(1:NROW(x), function(z) {min(x[z],y[z])})}
 rowMaxs <- function(x, y) {sapply(1:NROW(x), function(z) {max(x[z],y[z])})}

 hist(log(-minus), 80)
 hist(log(plus), 80)
 
 par(mfrow=c(2,1))
 hist(rowMins(log(-minus), log(plus)), 80, xlim=c(-4,10))
 hist(rowMaxs(log(-minus), log(plus)), 80, xlim=c(-4,10))

 sum(rowMins(log(-minus), log(plus))<0.5)/NROW(plus) # [1] 0.7165772
 sum(rowMaxs(log(-minus), log(plus))>0.5)/NROW(plus) # [1] 0.8914454

 all_data <- rbind(cbind(txn[log(plus)>0.5,], name="n", score=1, strand="+"), cbind(txn[log(-minus)>0.5,], name="n", score=1, strand="-"))
 
## Write it out and intersect clusters.
write.table(all_data, pipe(" sort-bed - > txn.strand.bed"), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

