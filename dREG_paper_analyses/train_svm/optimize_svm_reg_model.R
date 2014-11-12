## train_svm.R -- Trains an SVM to detect regulatory regions, and uses it to.
## scan the genome in 10bp windows.
##

setwd("/usr/projects/GROseq.parser/tss_detecter/")
require(featureDetector)

## Read PRO-seq data & positive regions.
ps_plus_path  <- "/usr/data/GROseq.parser/hg19/k562/proseq/K562_unt.subsamp10pct.bed.gz_plus.bw"  #"/usr/data/GROseq.parser/hg19/k562/proseq_celastrol_prelim/celastrol_proseq_0min_plus.bigWig"
ps_minus_path <- "/usr/data/GROseq.parser/hg19/k562/proseq/K562_unt.subsamp10pct.bed.gz_minus.bw" #"/usr/data/GROseq.parser/hg19/k562/proseq_celastrol_prelim/celastrol_proseq_0min_minus.bigWig"
GROcap_tss_bed <- read.table("/usr/projects/GROseq.parser/tss_new/hg19.k562.new_hmm2.bed", skip=1)
GROcap_merge_bed <- read.table("andrehmm.nostrand.merge.bed")


## Get informative positions for training.
inf_positions <- get_informative_positions(ps_plus_path, ps_minus_path, depth= 0, step=50, use_ANDOR=TRUE, use_OR=FALSE) ## Get informative positions.
NROW(inf_positions)

## Cut down here a little bit...
set.seed(61)
inf_positions <- inf_positions[sample(c(1:NROW(inf_positions)), 3000000, replace=FALSE), ]

#########################
## Make a test set.
n_eval <- 1000
tset <- get_test_set(inf_positions, GROcap_tss_bed, n_eval, allow= GROcap_merge_bed,  avoid_dist= 100)

## Search the space of possible representations.
gdm <- list()

## Test a set of reasonable window shapes.
gdm[[1]] <- genomic_data_model(window_sizes= c(10, 50, 500), half_nWindows= c(20, 20, 20)) # 0.953282388663967
gdm[[2]] <- genomic_data_model(window_sizes= c(10, 50, 500, 5000), half_nWindows= c(20, 20, 20, 20)) # 0.957300607287449
gdm[[3]] <- genomic_data_model(window_sizes= c(25, 50, 500, 5000), half_nWindows= c(20, 20, 20, 20)) # 0.955449392712551
gdm[[4]] <- genomic_data_model(window_sizes= c(10, 25, 50, 500, 5000), half_nWindows= c(20, 20, 20, 20, 20)) # 0.957096153846154 ## NO HELP.

gdm[[5]] <- genomic_data_model(window_sizes= c(10, 50, 500, 5000), half_nWindows= c(30, 20, 20, 20)) # 0.956994433198381 ## NO HELP.
gdm[[6]] <- genomic_data_model(window_sizes= c(10, 50, 500, 5000), half_nWindows= c(20, 30, 20, 20)) # 0.959086032388664 ## Slightly better (?!)
gdm[[7]] <- genomic_data_model(window_sizes= c(10, 50, 500, 5000), half_nWindows= c(20, 20, 30, 20)) # 0.957853238866397 ## NO HELP.
gdm[[8]] <- genomic_data_model(window_sizes= c(10, 50, 500, 5000), half_nWindows= c(20, 20, 20, 30)) # 0.956503036437247 ## NO HELP.

gdm[[9]] <- genomic_data_model(window_sizes= c(10, 50, 500, 5000), half_nWindows= c(10, 20, 20, 20)) # 0.957785425101214 ## Not much help?!  No hurt.
gdm[[10]] <- genomic_data_model(window_sizes= c(10, 50, 500, 5000), half_nWindows= c(20, 10, 20, 20)) # 0.956695344129555 ## HURT.
gdm[[11]] <- genomic_data_model(window_sizes= c(10, 50, 500, 5000), half_nWindows= c(20, 20, 10, 20)) # 0.955440789473684 ## HURT.
gdm[[12]] <- genomic_data_model(window_sizes= c(10, 50, 500, 5000), half_nWindows= c(20, 20, 20, 10)) # 0.955663461538462 ## HURT.

gdm[[13]] <- genomic_data_model(window_sizes= c(50, 500, 5000), half_nWindows= c(30, 20, 20)) # 0.953646761133603 ## MUCH worse.
gdm[[14]] <- genomic_data_model(window_sizes= c(10, 50, 500, 5000), half_nWindows= c(10, 30, 20, 20)) # 0.959619939271255 ## Best so far.
gdm[[15]] <- genomic_data_model(window_sizes= c(10, 50, 500, 5000), half_nWindows= c(15, 30, 20, 20)) # 0.958608299595142
gdm[[16]] <- genomic_data_model(window_sizes= c(10, 25, 50, 500, 5000), half_nWindows= c(10, 10, 30, 20, 20)) # 0.960782894736842 ## Best so far.

#pdf("test_GDM_sizes.pdf")
pdf("test_GDM_widths.trueLogistic.11-22-2013.pdf")
for(g in gdm) {
 print(g)
 set.seed(52)
 asvm <- regulatory_svm(g, ps_plus_path, ps_minus_path, inf_positions, GROcap_tss_bed, pdf_path= "roc_plot.and0.pdf", n_train=25000, n_eval=1000)
 save.image("Opt.current.svm.RData")
 
 ##############################################333
 ## Test model.
 pred_val<- eval_reg_svm(g, asvm, tset, ps_plus_path, ps_minus_path, batch_size= 10000)
 rocPlot <- roc.calc(GROcap_merge_bed, tset, pred_val)
 AUC <- roc.auc(rocPlot)
 roc.plot(rocPlot, main=paste("AUC=", AUC))

 print(paste("After training:", AUC))
 remove(asvm)
}
dev.off()
