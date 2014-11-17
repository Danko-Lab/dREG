## train_svm.R -- Trains an SVM to detect regulatory regions, and uses it to.
## scan the genome in 50bp windows.
##

#setwd("/usr/projects/GROseq.parser/tss_detecter/")
require(dREG)

## Read PRO-seq data.
ps_plus_path  <- "K562_unt.sort.bed.gz_plus.bw" #/usr/data/GROseq.parser/hg19/k562/proseq/
ps_minus_path <- "K562_unt.sort.bed.gz_minus.bw" #/usr/data/GROseq.parser/hg19/k562/proseq/

## Get positive regions.
GROcap_tss_bed <- read.table("hg19.k562.new_hmm2b.post2.dnase.bed") #"/usr/projects/GROseq.parser/tss_new/hg19.k562.new_hmm2.bed", skip=1)

## Train the SVM.
inf_positions <- get_informative_positions(ps_plus_path, ps_minus_path, depth= 0, step=50, use_ANDOR=TRUE, use_OR=FALSE) ## Get informative positions.
print(paste("Number of inf. positions: ", NROW(inf_positions)))

gdm <- genomic_data_model(window_sizes= c(10, 25, 50, 500, 5000), half_nWindows= c(10, 10, 30, 20, 20))
extra_enrich_bed <- read.table("GencodeMerge.IntersectOpStrand.bed")
allow_bed <- read.table("chromHmm.k562.enh.prom.bed.gz")
asvm <- regulatory_svm(gdm, ps_plus_path, ps_minus_path, inf_positions, GROcap_tss_bed, pdf_path= "roc_plot.and1.lgModel.pdf", n_train=50000, n_eval=5000, extra_enrich_bed= extra_enrich_bed, allow= allow_bed)

#remove(inf_positions, ps_plus_path, ps_minus_path, extra_enrich_bed, allow_bed)
save.image("asvm.getTrainSet.RData")

