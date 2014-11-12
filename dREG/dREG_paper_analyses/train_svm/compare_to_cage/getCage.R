require(bigWig)
require(featureDetector)
setwd("/usr/projects/GROseq.parser/tss_detecter/")

signal_threshold <- 8 ## Andre chose a threshold of 8 for both K562 and GM12878.

#pap_p <- "/usr/data/GROseq.parser/hg19/k562/cage_bp/riken_cage_nucleus_pap_plus.bigWig"
#pap_m <- "/usr/data/GROseq.parser/hg19/k562/cage_bp/riken_cage_nucleus_pap_minus.bigWig"
pap_p <- "/usr/data/GROseq.parser/hg19/gm12878/cage_bp/riken_cage_nucleus_pap_plus.bigWig"
pap_m <- "/usr/data/GROseq.parser/hg19/gm12878/cage_bp/riken_cage_nucleus_pap_minus.bigWig"

ip_p <- get_informative_positions(pap_p, depth=signal_threshold, window= 100, step=100)
ip_m <- get_informative_positions(pap_m, depth=signal_threshold, window= 100, step=100)

## Create a bed formatted dataset.
ip_p$name= rep("n")
ip_p$score= rep(0)
ip_p$strand= rep("+")

ip_m$name= rep("n")
ip_m$score= rep(0)
ip_m$strand= rep("-")

ip <- rbind(ip_p, ip_m)
ip$chromStart <- ip$chromStart -50
ip$chromEnds   <- ip$chromEnds +50

# Write a bed file.
options("scipen"=100, "digits"=4)
#write.table(ip, "cage_signal.k562.bed", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
write.table(ip, "cage_signal.gm12878.bed", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
