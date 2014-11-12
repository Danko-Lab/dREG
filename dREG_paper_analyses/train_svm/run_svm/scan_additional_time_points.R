## train_svm.R -- Trains an SVM to detect regulatory regions, and uses it to.
## scan the genome in 10bp windows.
##

setwd("/usr/projects/GROseq.parser/tss_detecter/")
require(featureDetector)

## Read PRO-seq data.
ps_plus_path  <- "/usr/data/GROseq.parser/hg19/k562/proseq_celastrol_prelim/celastrol_proseq_0min_plus.bigWig"
ps_minus_path <- "/usr/data/GROseq.parser/hg19/k562/proseq_celastrol_prelim/celastrol_proseq_0min_minus.bigWig"

## Get positive regions.
GROcap_tss_bed <- read.table("/usr/projects/GROseq.parser/tss_new/hg19.k562.new_hmm2.bed", skip=1)

## Train the SVM.
inf_positions <- get_informative_positions(ps_plus_path, ps_minus_path, depth= 4, use_OR=FALSE) ## Get informative positions.

gdm <- genomic_data_model(window_sizes= c(10, 50, 500, 5000), half_nWindows= c(20, 20, 20, 20))
asvm <- regulatory_svm(gdm, ps_plus_path, ps_minus_path, inf_positions, GROcap_tss_bed, pdf_path= "roc_plot.and4.pdf", n_train=15000, n_eval=2000)

save.image("RegulatorySVM.RData")

# setwd("/usr/projects/GROseq.parser/tss_detecter/")
# require(featureDetector)

# load("RegulatorySVM.RData")

## Get data for additional time points, and try to find regions.
ps_plus_path_20  <- "/usr/data/GROseq.parser/hg19/k562/proseq_celastrol_prelim/celastrol_proseq_20min_plus.bigWig"
ps_minus_path_20 <- "/usr/data/GROseq.parser/hg19/k562/proseq_celastrol_prelim/celastrol_proseq_20min_minus.bigWig"
inf_positions_20 <- get_informative_positions(ps_plus_path_20, ps_minus_path_20, depth= 4, use_OR=FALSE) ## Get informative positions.
pred_val_20<- eval_reg_svm(gdm, asvm, inf_positions_20, ps_plus_path_20, ps_minus_path_20)

final_data_20 <- data.frame(inf_positions_20, pred_val_20)
write.table(final_data_20, file="predictions_20m.bedGraph", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

ps_plus_path_80  <- "/usr/data/GROseq.parser/hg19/k562/proseq_celastrol_prelim/celastrol_proseq_80min_plus.bigWig"
ps_minus_path_80 <- "/usr/data/GROseq.parser/hg19/k562/proseq_celastrol_prelim/celastrol_proseq_80min_minus.bigWig"
inf_positions_80 <- get_informative_positions(ps_plus_path_80, ps_minus_path_80, depth= 4, use_OR=FALSE) ## Get informative positions.
pred_val_80<- eval_reg_svm(gdm, asvm, inf_positions_80, ps_plus_path_80, ps_minus_path_80)

final_data_80 <- data.frame(inf_positions_80, pred_val_80)
write.table(final_data_80, file="predictions_80m.bedGraph", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

## Get data for 0 and 5 min.
ps_plus_path_0  <- "/usr/data/GROseq.parser/hg19/k562/proseq_celastrol_prelim/celastrol_proseq_0min_plus.bigWig"
ps_minus_path_0 <- "/usr/data/GROseq.parser/hg19/k562/proseq_celastrol_prelim/celastrol_proseq_0min_minus.bigWig"
inf_positions_0 <- get_informative_positions(ps_plus_path_0, ps_minus_path_0, depth= 4, use_OR=FALSE) ## Get informative positions.
pred_val_0<- eval_reg_svm(gdm, asvm, inf_positions_0, ps_plus_path_0, ps_minus_path_0)

final_data_0 <- data.frame(inf_positions_0, pred_val_0)
write.table(final_data_0, file="predictions_0m.bedGraph", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

ps_plus_path_5  <- "/usr/data/GROseq.parser/hg19/k562/proseq_celastrol_prelim/celastrol_proseq_5min_plus.bigWig"
ps_minus_path_5 <- "/usr/data/GROseq.parser/hg19/k562/proseq_celastrol_prelim/celastrol_proseq_5min_minus.bigWig"
inf_positions_5 <- get_informative_positions(ps_plus_path_5, ps_minus_path_5, depth= 4, use_OR=FALSE) ## Get informative positions.
pred_val_5<- eval_reg_svm(gdm, asvm, inf_positions_5, ps_plus_path_5, ps_minus_path_5)

final_data_5 <- data.frame(inf_positions_5, pred_val_5)
write.table(final_data_5, file="predictions_5m.bedGraph", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

save.image("Predictions.RData")


###################################################################
## Test accuracy, and decide on a cutoff threshold.
require(featureDetector)

GROcap_tss_bed <- read.table("/usr/projects/GROseq.parser/tss_new/hg19.k562.new_hmm2.bed", skip=1)
GROcap_tss_bed[,2] <- GROcap_tss_bed[,2]-100 ## Anywhere within 100bp.
GROcap_tss_bed[,3] <- GROcap_tss_bed[,3]+100 ## Anywhere within 100bp.

pred <- read.table("predictions_0m.bedGraph.gz", skip=1)
indx <- sample(c(1:NROW(pred)), 50000, replace=FALSE)
roc_values <- roc.calc(GROcap_tss_bed[,c(1:3)], pred[indx,c(1:3)], pred$V4[indx])

roc.auc(roc_values)
pdf("Model.ROC.gt4.Celastrol-used.pdf")
 roc.plot(roc_values, main=paste("Celastrol; 0m Accuracy"))
 abline(v=0.05, col="red", lty="dotted")
dev.off()

#     FPR           TPR        THRESHOLD
#85 4.610041e-02 0.45440147  0.901184787
save_0m <- pred[pred$V4>0.9,]
write.table(save_0m, file="predictions_th0m.bedGraph", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

## Save other TPs.
pred_5 <- read.table("predictions_5m.bedGraph.gz", skip=1)
pred_5 <- pred_5[pred_5$V4>0.9,]
write.table(pred_5, file="predictions_th5m.bedGraph", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

pred_20 <- read.table("predictions_20m.bedGraph.gz", skip=1)
pred_20 <- pred_20[pred_20$V4>0.9,]
write.table(pred_20, file="predictions_th20m.bedGraph", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

pred_80 <- read.table("predictions_80m.bedGraph.gz", skip=1)
pred_80 <- pred_80[pred_80$V4>0.9,]
write.table(pred_80, file="predictions_th80m.bedGraph", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")


>>>>>>> 0cca60f32d18b3bc86e8c52fb7b0a1c602602066
