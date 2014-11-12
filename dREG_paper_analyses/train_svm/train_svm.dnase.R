## train_svm.R -- Trains an SVM to detect regulatory regions, and uses it to.
## scan the genome in 50bp windows.
##

#setwd("/usr/projects/GROseq.parser/tss_detecter/")
require(featureDetector)

## ALL DATA ORDER: K562, CD4, HeLa, MCF-7.

## Read PRO-seq data.
ps_plus_path  <- list("K562_unt.sort.bed.gz_plus.bw", "H-U_plus.bw", "HeLa.combined_plus.bw", "MCF7.unt.all_plus.bw")
ps_minus_path <- list("K562_unt.sort.bed.gz_minus.bw", "H-U_minus.bw", "HeLa.combined_minus.bw", "MCF7.unt.all_minus.bw")

## Get positive regions.
dnase <- lapply(c("wgEncodeUWDukeDnaseK562.fdr01peaks.hg19.bed", "dnase1.peaks_peaks.narrowPeak", "uw.merge.narrowPeak.bed", "wgEncodeAwgDnaseUwdukeMcf7UniPk.narrowPeak.gz"), function(x) {read.table(x)})
extra_enrich_bed <- lapply(c(1:4), function(x) {read.table("GencodeMerge.IntersectOpStrand.bed")})
allow_bed <- lapply(c("chromHmm.k562.enh.prom.bed.gz", "CD4.chromHMM.Ernst2010.hg19.Prom.Enh.bed", "helas3.chromhmm.prom.enh.ins.bed", "H3K4me1.liftOver.peaks.bed.gz"), function(x) {read.table(x)})

## Train the SVM.
inf_positions <- lapply(c(1:NROW(ps_plus_path)), function(x) {get_informative_positions(ps_plus_path[[x]], ps_minus_path[[x]], depth= 0, step=50, use_ANDOR=TRUE, use_OR=FALSE)}) ## Get informative positions.
print(paste("Number of inf. positions: ", NROW(inf_positions)))

gdm <- genomic_data_model(window_sizes= c(10, 25, 50, 500, 5000), half_nWindows= c(10, 10, 30, 20, 20))
asvm <- regulatory_svm(gdm, ps_plus_path, ps_minus_path, inf_positions, dnase, n_train=12500, n_eval=0, extra_enrich_bed= extra_enrich_bed, allow= allow_bed)

remove(inf_positions, ps_plus_path, ps_minus_path, extra_enrich_bed, allow_bed)
save.image("dnase1..asvm.RData")
