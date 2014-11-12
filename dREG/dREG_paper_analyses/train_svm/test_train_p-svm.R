## train_svm.R -- Trains an SVM to detect regulatory regions, and uses it to.
## scan the genome in 10bp windows.
##

#setwd("/usr/projects/GROseq.parser/tss_detecter/")
require(featureDetector)

## Read PRO-seq data.
ps_plus_path  <- "/usr/data/GROseq.parser/hg19/k562/proseq/K562_unt.sort.bed.gz_plus.bw" 
ps_minus_path <- "/usr/data/GROseq.parser/hg19/k562/proseq/K562_unt.sort.bed.gz_minus.bw" 

## Get positive regions.
GROcap_tss_bed <- read.table("/usr/projects/GROseq.parser/tss_new/hg19.k562.new_hmm2.bed", skip=1)

## Train the SVM.
inf_positions <- get_informative_positions(ps_plus_path, ps_minus_path, depth= 0, step=50, use_ANDOR=TRUE, use_OR=FALSE) ## Get informative positions.
print(paste("Number of inf. positions: ", NROW(inf_positions)))

gdm <- genomic_data_model(window_sizes= c(10, 25, 50, 500, 5000), half_nWindows= c(10, 10, 30, 20, 20)) 
extra_enrich_bed <- read.table("GencodeMerge.IntersectOpStrand.bed")
allow_bed <- read.table("chromHmm.k562.enh.prom.bed.gz")
psvm <- regulatory_svm(gdm, ps_plus_path, ps_minus_path, inf_positions, GROcap_tss_bed, svm_type="P_SVM", pdf_path=NULL, n_train=5000, n_eval=5000, extra_enrich_bed= extra_enrich_bed, allow= allow_bed)

save.image("p_svm.model_only.RData")

load("p_svm.model_only.RData")
pred_val<- eval_reg_svm(gdm, psvm, inf_positions[1:100,], ps_plus_path, ps_minus_path, batch_size= 10000, ncores=16)

final_data <- data.frame(inf_positions, pred_val)
options("scipen"=100, "digits"=4)
write.table(final_data, file="predictions.bedGraph", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

