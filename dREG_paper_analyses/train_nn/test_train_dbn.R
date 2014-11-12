## train_svm.R -- Trains an SVM to detect regulatory regions, and uses it to.
## scan the genome in 10bp windows.
##

## train_svm.R -- Trains an SVM to detect regulatory regions, and uses it to.
## scan the genome in 10bp windows.
##

setwd("/usr/projects/GROseq.parser/tss_detecter/")
require(featureDetector)

## Read PRO-seq data & positive regions.
## Read PRO-seq data.
ps_plus_path  <- "/usr/data/GROseq.parser/hg19/k562/proseq/K562_unt.sort.bed.gz_plus.bw"
ps_minus_path <- "/usr/data/GROseq.parser/hg19/k562/proseq/K562_unt.sort.bed.gz_minus.bw"
GROcap_tss_bed <- read.table("andre_hmm/hg19.k562.new_hmm2b.post2.bed")

## Get informative positions for training.
inf_positions <- get_informative_positions(ps_plus_path, ps_minus_path, depth= 0, step=50, use_ANDOR=TRUE, use_OR=FALSE) ## Get informative positions.
print(paste("Number of inf. positions: ", NROW(inf_positions)))

## Cut down here a little bit...
set.seed(61)
inf_positions <- inf_positions[sample(c(1:NROW(inf_positions)), 10000000, replace=FALSE), ]

## Search the space of possible representations.
gdm <- genomic_data_model(window_sizes= c(10, 25, 50, 500, 5000), half_nWindows= c(10, 10, 30, 20, 20))
extra_enrich_bed <- read.table("GencodeMerge.IntersectOpStrand.bed")
allow_bed <- read.table("chromHmm.k562.enh.prom.bed.gz")

# c(360,500,500,1000)
adbn <- dbn(layer_sizes= c(360,300,300,500), batch_size=100, cd_n=1, momentum_decay= 0.8, weight_cost= 1e-5, learning_rate=0.1)
adbn <- regulatory_dbn(gdm, adbn, ps_plus_path, ps_minus_path, inf_positions, GROcap_tss_bed, pdf_path= "dbn.roc.lgModel.pdf", n_train=100000, n_eval=5000, extra_enrich_bed= extra_enrich_bed, allow= allow_bed)

## Now scan all positions in the genome ...
#pred_val<- eval_reg_dbn(gdm, asvm, inf_positions, ps_plus_path, ps_minus_path, batch_size= 10000, ncores=16)
#
#final_data <- data.frame(inf_positions, pred_val)
#options("scipen"=100, "digits"=4)
#write.table(final_data, file="predictions.bedGraph", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")


save.image("adbn.1.RData")
