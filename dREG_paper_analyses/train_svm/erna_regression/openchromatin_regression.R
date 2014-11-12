##
##
## Goal: Regress the presence/ absence of eRNA against 
##       a set of transcription factors.  Identify the 
##       set of TFs that drive eRNA expression, and those
##       that do not.
##

fold_cv <- 0.9
dreg_max_neg<- 0.5   ## MAX of the negative samples.
#fold_cv <- 0.9

## Read in DNAse-1 peaks
## V7 is signal ... see here: http://genome.ucsc.edu/FAQ/FAQformat.html#format12
#dnase1 <- read.table("/usr/data/GROseq.parser/hg19/k562/dnase/wgEncodeOpenChromDnaseK562PkV2.narrowPeak.gz")

## dREG
#dreg <- read.table("k562.scores.bed")

peaks <- read.table("overlaps.dnase1.chromhmm.tsv")
tfs <- read.table("/usr/data/GROseq.parser/hg19/k562/tf_peaks/tf_names_and_files.txt")
colnames(peaks) <- c("chrom", "start", "end", "DNAse-1.peakscore", "dreg", as.character(tfs$V1))

## Fraction of sites in regions without DNAse-1
colSums(peaks[peaks$"DNAse-1.peakscore"=="NAN",c(6:NCOL(peaks))])/colSums(peaks[,c(6:NCOL(peaks))]) ## Without DNAse-1.
wdnase <- colSums(peaks[!peaks$"DNAse-1.peakscore"=="NAN",c(6:NCOL(peaks))])/colSums(peaks[,c(6:NCOL(peaks))]) ## With DNAse-1 peak.
require(boot)
boot_wd <- boot(data= peaks, R= 1000, statistic= function(a, i){
  colSums(a[i,][!(a$"DNAse-1.peakscore"=="NAN")[i],c(6:NCOL(peaks))])/colSums(a[i,][,c(6:NCOL(peaks))])
})
error_wd <- sapply(1:NROW(boot_wd$t0) , function(x) {sd(boot_wd$t[,x], na.rm=TRUE)})

## Divide into train and test sets.
use <- is.na(peaks$dreg) | peaks$dreg<dreg_max_neg #rep(TRUE, NROW(peaks)) #peaks$dreg>dreg_th | peaks$dreg<dreg_max_neg | is.na(peaks$dreg) ## Remove peaks just below threshold.

train <- sample(c(1:NROW(peaks))[use], sum(use)*fold_cv)
test  <- rep(TRUE, NROW(peaks)); test[train] <- FALSE; test[!use] <- FALSE; test <- which(test)

## Now create a regression.
data_df <- cbind(data.frame(y= !(peaks$"DNAse-1.peakscore" == "NAN"), peaks[,c(6:NCOL(peaks))]))
fm <- glm(y~., family=binomial, data=data_df[train,])			## Full model with TFs.

## How well does the model explain the variation?  How accurately can we predict which DNAse-1 sites will have eRNAs?
scores_fm <- predict(fm, data_df[test,]) 

require(featureDetector)
roc_fm <- logreg.roc.calc(data_df$y[test], scores_fm)
roc.auc(roc_fm)

pdf("chromhmmdnase.pdf")
 source("~/src/featureDetector/test_functions/figures/histplot.R")
 cd.barplot(wdnase, error_wd, names(wdnase), "black")
 
 roc.plot(roc_fm, xlim=c(0,1), ylim=c(0,1), col="dark red")

## Now use bootstrap to test the regression coefficients difference from 0.
require(boot)
bb <- boot(data= data_df[use,], R= 1000, statistic= function(a, i) {
  glm(y~., family=binomial, data=a[i,])$coefficients
})
save.image("~/src/featureDetector/test_functions/train_svm/erna_regression/opensite_regression.RData")

std.error <- sapply(1:NROW(bb$t0) , function(x) {sd(bb$t[,x], na.rm=TRUE)})
sig <- sapply(1:NROW(bb$t0) , function(x) {!xor(quantile(bb$t[,x], 0.025, na.rm=TRUE)>0, quantile(bb$t[,x], 0.975, na.rm=TRUE)>0)})
draw <- !is.na(bb$t0)

#minLt <- sapply(1:NROW(bb$t0) , function(x) {min(bb$t[,x])<-1e10})

source("~/src/featureDetector/test_functions/train_svm/erna_regression/erna_drawbars.R")
drawBars(bb$t0[draw], std.error[draw], names(bb$t0)[draw])
drawBars(bb$t0[sig&draw], std.error[sig&draw], names(bb$t0)[sig&draw])

dev.off()

