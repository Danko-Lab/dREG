## train_svm.R -- Trains an SVM to detect regulatory regions, and uses it to.
## scan the genome in 10bp windows.
##

#setwd("/usr/projects/GROseq.parser/tss_detecter/")
require(featureDetector)

## Read PRO-seq data.
ps_plus_path  <- "K562_unt.sort.bed.gz_plus.bw" #/usr/data/GROseq.parser/hg19/k562/proseq/
ps_minus_path <- "K562_unt.sort.bed.gz_minus.bw" #/usr/data/GROseq.parser/hg19/k562/proseq/

## Get positive regions.
GROcap_tss_bed <- read.table("hg19.k562.new_hmm2b.post2.bed") #"/usr/projects/GROseq.parser/tss_new/hg19.k562.new_hmm2.bed", skip=1)

## Train the SVM.
inf_positions <- get_informative_positions(ps_plus_path, ps_minus_path, depth= 0, step=50, use_ANDOR=TRUE, use_OR=FALSE) ## Get informative positions.
print(paste("Number of inf. positions: ", NROW(inf_positions)))

gdm <- genomic_data_model(window_sizes= c(10, 25, 50, 500, 5000), half_nWindows= c(10, 10, 30, 20, 20)) 
extra_enrich_bed <- read.table("GencodeMerge.IntersectOpStrand.bed")
allow_bed <- read.table("chromHmm.k562.enh.prom.bed.gz")
asvm <- regulatory_svm(gdm, ps_plus_path, ps_minus_path, inf_positions, GROcap_tss_bed, pdf_path= "roc_plot.and1.lgModel.pdf", n_train=50000, n_eval=5000, extra_enrich_bed= extra_enrich_bed, allow= allow_bed)

save.image("asvm.model_only.RData")

## Now scan all positions in the genome ...
inf_positions <- get_informative_positions(ps_plus_path, ps_minus_path, depth= 0, step=50, use_ANDOR=TRUE, use_OR=FALSE) ## Get informative positions.

pred_val<- eval_reg_svm(gdm, asvm, inf_positions, ps_plus_path, ps_minus_path, batch_size= 10000, ncores=16)

final_data <- data.frame(inf_positions, pred_val)
options("scipen"=100, "digits"=4)
write.table(final_data, file="predictions.bedGraph", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

save.image("asvm.RData")

q("no")
### FORCE QUIT ####################

###################################
## Analyze the trained model ... 

##  Finally, we'll merge consecutive positive values and intersect with Andre's bigWig.
final_data <- read.table("k562.bed.gz")

require(featureDetector)
GROcap_tss_bed <- read.table("andre_hmm/hg19.k562.new_hmm2b.post2.bed") #read.table("/usr/projects/GROseq.parser/tss_new/hg19.k562.new_hmm2.bed", skip=1)

dnase1 <- read.table("/usr/data/GROseq.parser/hg19/k562/dnase/wgEncodeOpenChromDnaseK562PkV2.narrowPeak.gz")

#zcat /usr/data/GROseq.parser/hg19/k562/chromhmm/wgEncodeBroadHmmK562HMM.bed.gz | grep "Promoter\|Enhancer" | sort-bed - | bedops --merge - > /usr/projects/GROseq.parser/tss_detecter/k562.chromhmm.promoter.enhancer.merge.bed
chromhmm <- read.table("k562.chromhmm.promoter.enhancer.merge.bed")

rocPlot <- roc.calc(GROcap_tss_bed, final_data, final_data[,4], filterPossible=TRUE)

roc.auc(rocPlot)
roc.plot(rocPlot)

## 
pdf("MetaPlots.pdf")

## Histograms of both ...
# AFTER filtering...
par(mfrow=c(4,1))

hist(final_data$V3-final_data$V2, seq(0,10000,50), xlim=c(0,3000), main="TSS Detector", xlab="TSS size [bp]")
hist(GROcap_tss_bed$V3-GROcap_tss_bed$V2, seq(0,10000,50), xlim=c(0,3000), main="GRO-cap", xlab="TSS size [bp]")
hist(dnase1$V3-dnase1$V2, seq(0,1000000,50), xlim=c(0,3000), main="DNAse-1", xlab="TSS size [bp]")
hist(chromhmm$V3-chromhmm$V2, seq(0,1000000,200), xlim=c(0,3000), main="chromHMM", xlab="TSS size [bp]")

par(mfrow=c(1,2))

## Make meta-plots of ...
halfWindow=5000
step=100


par(mfrow=c(2,1))

## PRO-seq.
require(bigWig)
proseq_plus_bw <- load.bigWig("/usr/data/GROseq.parser/hg19/k562/proseq/K562_unt.sort.bed.gz_plus.bw")
proseq_minus_bw <- load.bigWig("/usr/data/GROseq.parser/hg19/k562/proseq/K562_unt.sort.bed.gz_minus.bw")
proseq_plus_meta<- meta.subsample(bed= final_data, bigWig.plus= proseq_plus_bw, halfWindow=halfWindow, step=step)
proseq_minus_meta<- meta.subsample(bed= final_data, bigWig.plus= proseq_minus_bw, halfWindow=halfWindow, step=step)
meta.plot.GROseq(proseq_plus_meta, proseq_minus_meta, step=step, main="PRO-seq")

proseq_AM_meta <- meta.subsample(bed= GROcap_tss_bed, bigWig.plus= proseq_plus_bw, bigWig.minus= proseq_minus_bw, halfWindow=halfWindow, step=step)
meta.plot.GROseq(proseq_AM_meta, proseq_AM_meta, step=step, main="PRO-seq ... Andre's peaks.")

## GRO-cap.
grocap_plus_bw <- load.bigWig("/usr/data/GROseq.parser/hg19/k562/groseq_tss/groseq_tss_wTAP_plus.bigWig")
grocap_minus_bw <- load.bigWig("/usr/data/GROseq.parser/hg19/k562/groseq_tss/groseq_tss_wTAP_minus.bigWig")
grocap_plus_meta<- meta.subsample(bed= final_data, bigWig.plus= grocap_plus_bw, halfWindow=halfWindow, step=step)
grocap_minus_meta<- meta.subsample(bed= final_data, bigWig.plus= grocap_minus_bw, halfWindow=halfWindow, step=step)
meta.plot.GROseq(grocap_plus_meta, grocap_minus_meta, step=step, main="GRO-cap")

grocap_AM_meta <- meta.subsample(bed= GROcap_tss_bed, bigWig.plus= grocap_plus_bw, bigWig.minus= grocap_minus_bw, halfWindow=halfWindow, step=step)
meta.plot.GROseq(grocap_AM_meta, grocap_AM_meta, step=step, main="GRO-cap ... Andre's peaks.")

par(mfrow=c(2,2))

## Histone marks.
H3K4me1_bw <- load.bigWig("/usr/data/GROseq.parser/hg19/k562/histones/wgEncodeBroadHistoneK562H3k4me1StdSig.bigWig")
H3K4me1_meta<- meta.subsample(bed= final_data, bigWig.plus= H3K4me1_bw, halfWindow=halfWindow, step=step)
meta.plot(H3K4me1_meta, step=step, main="H3K4me1")

H3K4me3_bw <- load.bigWig("/usr/data/GROseq.parser/hg19/k562/histones/wgEncodeBroadHistoneK562H3k4me3StdSig.bigWig")
H3K4me3_meta<- meta.subsample(bed= final_data, bigWig.plus= H3K4me3_bw, halfWindow=halfWindow, step=step)
meta.plot(H3K4me3_meta, step=step, main="H3K4me3")

H3K9ac_bw <- load.bigWig("/usr/data/GROseq.parser/hg19/k562/histones/wgEncodeBroadHistoneK562H3k9acStdSig.bigWig")
H3K9ac_meta<- meta.subsample(bed= final_data, bigWig.plus= H3K9ac_bw, halfWindow=halfWindow, step=step)
meta.plot(H3K9ac_meta, step=step, main="H3K9ac")

H3K27ac_bw <- load.bigWig("/usr/data/GROseq.parser/hg19/k562/histones/wgEncodeBroadHistoneK562H3k27acStdSig.bigWig")
H3K27ac_meta<- meta.subsample(bed= final_data, bigWig.plus= H3K27ac_bw, halfWindow=halfWindow, step=step)
meta.plot(H3K27ac_meta, step=step, main="H3K27ac")

dev.off()

## Select non-overlapping elements
#cat /usr/projects/GROseq.parser/tss_new/hg19.k562.new_hmm2.bed | grep "^t" -v > tmp.andre.tss.tmp.bed
#overlapSelect -nonOverlapping tmp.andre.tss.tmp.bed predictions.bed non-tss.bed

## Make meta-plots of histone marks ...
