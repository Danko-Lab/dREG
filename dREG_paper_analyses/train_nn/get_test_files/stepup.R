require(bigWig)

gcplus  <- load.bigWig("/usr/data/GROseq.parser/hg19/k562/groseq_tss/groseq_tss_wTAP_plus.bigWig")
gcminus <- load.bigWig("//usr/data/GROseq.parser/hg19/k562/groseq_tss/groseq_tss_wTAP_minus.bigWig")

grocap <- read.table("/usr/projects/GROseq.parser/tss_detecter/andre_hmm/hg19.k562.new_hmm2b.post2.bed")

## Get the expression of GRO-cap in each cluster.  
grocap[,5] <- bedQuery.bigWig(grocap, gcplus, gcminus)

## Intersect with the 5' end of known genes.

## Write it out and intersect clusters.
write.table(grocap, "grocap_counts.")

