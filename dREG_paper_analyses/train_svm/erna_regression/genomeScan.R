## Tests dREG as a 
##
require(rtfbsdb)
surroundSize <- 200

## Feed rtfbsdb the entire genome.
offset_dist <- 250
#system("twoBitInfo /gbdb/hg19/hg19.2bit hg19.chromInfo")
chromInfo <- read.table("~/storage/data/hg19/hg19.chromInfo")
chromInfo <- chromInfo[grep("_|chrM|chrY|chrX", chromInfo[,1], invert=TRUE),] ## Remove random, M, x, and y
chromInfo <- data.frame(chrom=chromInfo[,1], chromStart=rep(0)+offset_dist, chromEnd=(chromInfo[,2]-1-offset_dist))

## Now get it going!
## Initail site scane completed using PIQ.
#NRF1.sites <- scan_rtfbs("NRF1", chromInfo, 
#				"/usr/data/GROseq.parser/pwm_data/jolma/teal/NRF1.YGCGCATGCGCN.pwm", 
#				twoBit_path= "/gbdb/hg19/hg19.2bit", 
#				return_posteriors=FALSE, 
#				threshold = 2)
#NROW(NRF1.sites)
#write.table(NRF1.sites, pipe("sort-bed - | starch - > NRF1.starch"), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
#
#ELF1.sites <- scan_rtfbs("ELF1", chromInfo, 
#				"/usr/data/GROseq.parser/pwm_data/jolma/teal/ELF1.AACCCGGAAGTR.pwm", 
#				twoBit_path= "/gbdb/hg19/hg19.2bit", 
#				return_posteriors=FALSE, 
#				threshold = 2)
#NROW(ELF1.sites)
#write.table(ELF1.sites, pipe("sort-bed - | starch - > ELF1.starch"), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
#
#SP1.sites <- scan_rtfbs("SP1", chromInfo, 
#				"/usr/data/GROseq.parser/pwm_data/jolma/teal/SP1.GCCMCGCCCMC.pwm", 
#				twoBit_path= "/gbdb/hg19/hg19.2bit", 
#				return_posteriors=FALSE, 
#				threshold = 2)
#NROW(SP1.sites)
#write.table(SP1.sites, pipe("sort-bed - | starch - > SP1.starch"), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
#
#MAX.sites <- scan_rtfbs("MAX", chromInfo, 
#				"/usr/data/GROseq.parser/pwm_data/jolma/teal/MAX.NNCACGTGNN.pwm", 
#				twoBit_path= "/gbdb/hg19/hg19.2bit", 
#				return_posteriors=FALSE, 
#				threshold = 2)
#NROW(MAX.sites)
#write.table(MAX.sites, pipe("sort-bed - | starch - > MAX.starch"), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

## Parse binding sites/ peaks w/ bedops.
system("bash genomeScan.bsh")

analyzeData <- function(filename) {
  ## Read data.
  test_tf <- read.table(filename); test_tf$V12[is.na(test_tf$V12)] <- 0; test_tf <- test_tf[grep("_|chrM|chrY|chrX", test_tf$V1, invert=TRUE),]
  left <- test_tf[,2]-surroundSize; left[left < 0] <- 0 ## Should fix bigWig error trying to cout coordinates <0

  ## Get DNAse-1/ GRO-seq counts in the same region.
  require(bigWig)
  dnase <- load.bigWig("/home/cgd24/storage/data/hg19/k562/dnase/wgEncodeUwDgfK562Sig.bigWig")
  gsPlus <- load.bigWig("/home/cgd24/storage/data/hg19/k562/proseq/K562_unt.sort.bed.gz_plus.bw")
  gsMinus <- load.bigWig("/home/cgd24/storage/data/hg19/k562/proseq/K562_unt.sort.bed.gz_minus.bw")
  test_tf <- cbind(test_tf, 
					DNAse= bed.region.bpQuery.bigWig(dnase, data.frame(test_tf[,1], left, test_tf[,3]+surroundSize)), 
					GROseq_fwd= bed.region.bpQuery.bigWig(gsPlus, data.frame(test_tf[,1], left, test_tf[,3]+surroundSize)),
					GROseq_rev= abs(bed.region.bpQuery.bigWig(gsMinus, data.frame(test_tf[,1], left, test_tf[,3]+surroundSize))))
  #cor.test(test_tf$V7, test_tf$V8, method="spearman"); cor.test(test_tf$DNAse, test_tf$V8, method="spearman"); cor.test(test_tf$GROseq_fwd+test_tf$GROseq_rev, test_tf$V8, method="spearman")

  ## Create a linear model of binding.
  train <- sample(c(1:NROW(test_tf)), NROW(test_tf)*0.5)
  test  <- rep(TRUE, NROW(test_tf)); test[train] <- FALSE; test <- which(test)

  test_tf.lm <- glm(V13~V12+DNAse+GROseq_fwd+GROseq_rev+V8, family=binomial, data=test_tf[train,]); test_tf <- cbind(test_tf, test_tf.glm= predict(test_tf.lm, test_tf))
  test_gs.lm <- glm(V13~V12+GROseq_fwd+GROseq_rev, family=binomial, data=test_tf[train,]); test_tf <- cbind(test_tf, test_gs.glm= predict(test_gs.lm, test_tf))

  ## TRUNCATE TO TEST SET
  test_tf_eval <- test_tf[test,]
	
  ## Create ROC plot comparison.
  require(featureDetector)

  test_tf.roc <- logreg.roc.calc(test_tf_eval$V13, test_tf_eval$V12) #DREG
  test_tfD.roc <- logreg.roc.calc(test_tf_eval$V13, test_tf_eval$DNAse) #DNase (simple).
  test_tfPIQ.roc <- logreg.roc.calc(test_tf_eval$V13, test_tf_eval$V11) #PIQ PPV guestimate.
  test_tfM.roc <- logreg.roc.calc(test_tf_eval$V13, test_tf_eval$V8) # Motif
  test_tfGLM.roc <- logreg.roc.calc(test_tf_eval$V13, test_tf_eval$test_tf.glm) 
  test_tfGSGLM.roc <- logreg.roc.calc(test_tf_eval$V13, test_tf_eval$test_gs.glm)

  print(paste("GRO-seq: ", roc.auc(test_tf.roc)))
  print(paste("DNAse-1: ", roc.auc(test_tfD.roc)))
  print(paste("PIQ: ", roc.auc(test_tfPIQ.roc)))
  print(paste("Motif: ", roc.auc(test_tfM.roc)))
  print(paste("GLM: ", roc.auc(test_tfGLM.roc)))
  print(paste("GRO-seq GLM: ", roc.auc(test_tfGSGLM.roc)))
  
  print(paste("GRO-seq Performance: ", (roc.auc(test_tf.roc)-roc.auc(test_tfPIQ.roc))/roc.auc(test_tfPIQ.roc),"(PIQ)", (roc.auc(test_tf.roc)-roc.auc(test_tfD.roc))/roc.auc(test_tfD.roc),"(DNAse)", (roc.auc(test_tf.roc)-roc.auc(test_tfM.roc))/roc.auc(test_tfM.roc), "(Motif)"))
	
  roc.plot(test_tf.roc, col="#e30000", xlim=c(0,1), ylim=c(0,1), main=filename)
  par(new=TRUE); roc.plot(test_tfD.roc, col="black", xlim=c(0,1), ylim=c(0,1))
  par(new=TRUE); roc.plot(test_tfPIQ.roc, col="dark green", xlim=c(0,1), ylim=c(0,1))
  par(new=TRUE); roc.plot(test_tfM.roc, col="gray", xlim=c(0,1), ylim=c(0,1))
  par(new=TRUE); roc.plot(test_tfGLM.roc, col="red", xlim=c(0,1), ylim=c(0,1), lwd=2)
  
  return(test_tf)
}

pdf("TFROCS.pdf")
 seqLogo(exp(t(read.motif("/home/cgd24/storage/data/pwm_data/jolma/teal/NRF1.YGCGCATGCGCN.pwm", header=TRUE))), newpage=TRUE)
 nrf1_ <- analyzeData("NRF1.scores.bed")
 
 seqLogo(exp(t(read.motif("/home/cgd24/storage/data/pwm_data/jolma/teal/ELF1.AACCCGGAAGTR.pwm", header=TRUE))), newpage=TRUE)
 elf1_ <- analyzeData("ELF1.scores.bed")

 seqLogo(exp(t(read.motif("/home/cgd24/storage/data/pwm_data/jolma/teal/SP1.GCCMCGCCCMC.pwm", header=TRUE))), newpage=TRUE)
 sp1_ <- analyzeData("SP1.scores.bed")

 seqLogo(exp(t(read.motif("/home/cgd24/storage/data/pwm_data/jolma/teal/MAX.NNCACGTGNN.pwm", header=TRUE))), newpage=TRUE)
 max_ <- analyzeData("MAX.scores.bed")

 seqLogo(exp(t(read.motif("/home/cgd24/storage/data/pwm_data/jolma/teal/SP1.GCCMCGCCCMC.pwm", header=TRUE))), newpage=TRUE)
 sp1jolma_ <- analyzeData("SP1JASPAR.scores.bed")
dev.off()
