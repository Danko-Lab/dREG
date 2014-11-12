##
## train_nnet.R -- trains a traditional neural network using the 
## backpropagation algorithm.
##

setwd("/usr/projects/GROseq.parser/tss_detecter/")

require(bigWig)
require(nnet)
require(featureDetector)
source("roc.calc.R")

###################################
## Vars. 
n.examp <- 10000
step.size   <- 20

###################################
## Compile data for positive/ neg.

## Train on whole genome, using Andre's GRO-cap TSS 'predictions' as the true set.
GROcap_tss_bed <- read.table("/usr/projects/GROseq.parser/tss_new/hg19.k562.new_hmm2.bed", skip=1)
rand_train_tss_bed <- GROcap_tss_bed[sort(sample(c(1:NROW(GROcap_tss_bed)), n.examp, replace=FALSE)),c(1:3)]  ## Train from just he plus strand ... we have enough, and the negative strand expected to look the same.

## Get transcription units.
Gencode <- read.table("/usr/projects/GROseq.parser/annotations/gencode.comprehensive.bed", header=FALSE, skip=1)
Gencode <- Gencode[Gencode[,11]== "protein_coding",]
rand_train_notss_bed <- Gencode[sort(sample(c(1:NROW(Gencode)), n.examp, replace=FALSE)),c(1:3)]

# Combine.
rand_train_bed <- rbind(rand_train_tss_bed, rand_train_notss_bed)
rand_train_classes <- c(rep(1, NROW(rand_train_tss_bed)), rep(0, NROW(rand_train_notss_bed)))

# Now make a large test set for evaluating preformance. 
# IDEA: Try to replace this with manually curated set.  Optimize w.r.t. the manually curated set.
rand_test_tss_bed <- GROcap_tss_bed[sort(sample(c(1:NROW(Gencode)), n.examp, replace=FALSE)),c(1:3)]  ## Evaluate.
rand_test_notss_bed<- Gencode[sort(sample(c(1:NROW(GROcap_tss_bed)), n.examp, replace=FALSE)),c(1:3)]
rand_test_bed <- rbind(rand_test_tss_bed, rand_test_notss_bed)
rand_test_classes <- c(rep(1, NROW(rand_test_tss_bed)), rep(0, NROW(rand_test_notss_bed)))

###################################
## Functions to train and evaluate.

## Optimizes number of hidden nodes, and returns the optimal model.
opt.nnet <- function(gs_plus, gs_minus, x_train_bed, x_predict_bed, y_train, y_predict= y_train, nhidden= c(5, 50, 100, 200), ...) {

 print("Collecting training data.")
 x_train <- read_genomic_data(x_train_bed, gs_plus, gs_minus)
 
 print("Collecting predicted data.")
 x_predict <- read_genomic_data(x_predict_bed, gs_plus, gs_minus)

 print("Training a model")
# Fit a feed-forard neural networks w/ each hidden node size.
 mod <- NULL
 acc <- -1
 for(i in nhidden) {
  n <- nnet(y= y_train, x= x_train, size=i, MaxNWts= 500000, linout= TRUE) ## size -> num. nodes in the hidden layer.
  pred <- predict(n, x_predict)
  pred[pred >  0.5] <- 1
  pred[pred <= 0.5] <- 0
  pred_acc <- (sum(pred == y_predict)/NROW(y_predict))
  if(acc < pred_acc) mod <- n
  print(paste("#:", i, "A:", pred_acc))
 }
 
 ## Plot a ROC plot.
 roc_values <- logreg.roc.calc(y_predict, predict(mod, x_predict))
 roc.plot(roc_values, ...)

# Return the best performing neural network.
 return(mod)
}

require(e1071)
opt.svm <- function(gs_plus, gs_minus, x_train_bed, x_predict_bed, y_train, y_predict= y_train, print_raw_data= TRUE, ...) {

 print("Collecting training data.")
 x_train <- read_genomic_data(x_train_bed, gs_plus, gs_minus)
 
 print("Collecting predicted data.")
 x_predict <- read_genomic_data(x_predict_bed, gs_plus, gs_minus)

 asvm <- svm(x_train, y_train)
 pred <- predict(asvm, x_predict)
 pred[pred >  0.5] <- 1
 pred[pred <= 0.5] <- 0
 pred_acc <- (sum(pred == y_predict)/NROW(y_predict))
 print(paste("A:", pred_acc))

 roc_values <- logreg.roc.calc(y_predict, predict(asvm, x_predict))
 roc.plot(roc_values, ...)
 
 if(print_raw_data) {
  plot(colSums(x_train[y_train == 1,]), ylab="Training data", type="l", ...)
  points(colSums(x_train[y_train == 0,]), col="gray", type="l", ...)
  plot(colSums(x_predict[y_predict == 1,]), ylab="Prediction data", type="l", ...)
  points(colSums(x_predict[y_predict == 0,]), col="gray", type="l", ...)
 }

 return(asvm)
}


###################################
## Data source information.   Trains a model for each data type.
pdf("GROseq.plots.pdf")
## TSS. 
## GRO-cap.
groCap.plus.path  <- "/usr/data/GROseq.parser/hg19/k562/groseq_tss/groseq_tss_wTAP_plus.bigWig" 
groCap.minus.path <- "/usr/data/GROseq.parser/hg19/k562/groseq_tss/groseq_tss_wTAP_minus.bigWig"

#groCap.minus.model <- opt.nnet(datapath= groCap.minus.path, rand.train.bed, rand.test.bed, rand.train.classes, main="GRO-cap Minus")

## Cage.
#cage.plus.path  <- "/usr/data/GROseq.parser/hg19/k562/cage/wgEncodeRikenCageK562NucleusPamPlusSignal.bigWig"
#cage.minus.path <- "/usr/data/GROseq.parser/hg19/k562/cage/wgEncodeRikenCageK562NucleusPamMinusSignal.bigWig"

## GRO-seq
gs_plus  <- "/usr/data/GROseq.parser/hg19/k562/groseq/groseq_plus.bigWig" 
gs_minus <- "/usr/data/GROseq.parser/hg19/k562/groseq/groseq_minus.bigWig" 
proseq_svm_model <- opt.svm(gs_plus= gs_plus, gs_minus= gs_minus, x_train_bed= rand_train_bed, x_predict_bed= rand_test_bed, y_train= rand_train_classes, main="GRO-seq SVM")
proSeq_model <- opt.nnet(gs_plus= gs_plus, gs_minus= gs_minus, x_train_bed= rand_train_bed, x_predict_bed= rand_test_bed, y_train= rand_train_classes, main="GRO-seq Nerual Net")

## PRO-seq.
ps_plus  <- "/usr/data/GROseq.parser/hg19/k562/proseq_celastrol_prelim/celastrol_proseq_0min_plus.bigWig"
ps_minus <- "/usr/data/GROseq.parser/hg19/k562/proseq_celastrol_prelim/celastrol_proseq_0min_minus.bigWig"
proseq_svm_model <- opt.svm(gs_plus= ps_plus, gs_minus= ps_minus, x_train_bed= rand_train_bed, x_predict_bed= rand_test_bed, y_train= rand_train_classes, main="PRO-seq SVM")
proSeq_model <- opt.nnet(gs_plus= ps_plus, gs_minus= ps_minus, x_train_bed= rand_train_bed, x_predict_bed= rand_test_bed, y_train= rand_train_classes, main="PRO-seq Nerual Net")

dev.off()

# QUIT ##############################################################################################3
q("no")
# QUIT ##############################################################################################3
# QUIT ##############################################################################################3
# QUIT ##############################################################################################3
# QUIT ##############################################################################################3

## Histone mods.
histone.dir <- "/usr/data/GROseq.parser/hg19/k562/histones/"
histone.files <- paste(histone.dir, c("wgEncodeBroadHistoneK562H3k4me3StdSig.bigWig", "wgEncodeBroadHistoneK562H3k4me1StdSig.bigWig",
					"wgEncodeBroadHistoneK562H3k9acStdSig.bigWig", "wgEncodeBroadHistoneK562H3k27acStdSig.bigWig"), sep="")#	"", "")  #bw.files[["h3k4me3"]] <- load.bigWig(histone.files[1])#histone.h3k4me3 <- #bw.files[["h3k4me1"]] <- load.bigWig(histone.files[2])#histone.h3k4me1 <- #bw.files[["h3k9ac"]]  <- load.bigWig(histone.files[3])#histone.H3k9ac  <- #bw.files[["h3k27ac"]] <- load.bigWig(histone.files[4])#histone.H3k27ac <- 					

histone.h3k4me3.model <- opt.nnet(datapath= histone.files[1], rand.train.bed, rand.test.bed, rand.train.classes, main="H3K4me3")
histone.H3k4me1.model <- opt.nnet(datapath= histone.files[2], rand.train.bed, rand.test.bed, rand.train.classes, main="H3K4me1")
histone.H3k9ac.model  <- opt.nnet(datapath= histone.files[3], rand.train.bed, rand.test.bed, rand.train.classes, main="H3K9ac")
histone.H3k27ac.model <- opt.nnet(datapath= histone.files[4], rand.train.bed, rand.test.bed, rand.train.classes, main="H3K27ac")

## DNAse-1.
dnase1.path <- "/usr/data/GROseq.parser/hg19/k562/dnase/wgEncodeUwDgfK562Sig.bigWig"
dnase1.model <- opt.nnet(datapath= dnase1.path, rand.train.bed, rand.test.bed, rand.train.classes, main="DNAse-1")

#save.image("tss.models.RData")

###################################
## Combine nnet models and evaluates performance ...

## Nerual network, of basepair resolution data in surrounding region.
## Learn from: 
#require(nnet) ## Fit single-hidden-layer neural network, possibly with skip-layer connections.
#n <- nnet(y= classes.bin, x= rbind(bw.pos.data[[1]], bw.neg.data[[1]]), size=100, MaxNWts= 15000) ## size -> num. nodes in the hidden layer.
#predict(n, rbind(bw.pos.data[[1]], bw.neg.data[[1]]))

