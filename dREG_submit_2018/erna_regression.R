##
##
## Goal: Regress the presence/ absence of eRNA against
##       a set of transcription factors.  Identify the
##       set of TFs that drive eRNA expression, and those
##       that do not.
##

args <- commandArgs(trailingOnly=TRUE)
pre   <- args[1]

fold_cv <- 0.9
sl <- (1e-3)/2


library(parallel)

file.G1 <- "../new-rf-201803/G1/G1.dREG.peak.score.bed.gz"
file.G2 <- "../new-rf-201803/G2/G2.dREG.peak.score.bed.gz"
file.G3 <- "../new-rf-201803/G3/G3.dREG.peak.score.bed.gz"
file.G4 <- "../new-rf-201803/G4/G4.dREG.peak.score.bed.gz"
file.G5 <- "../new-rf-201803/G5/G5.dREG.peak.score.bed.gz"
file.G6 <- "../new-rf-201803/G6/G6.dREG.peak.score.bed.gz"
file.G7 <- "../new-rf-201803/G7/G7.dREG.peak.score.bed.gz"

file.TF.chipseq <- "/fs/cbsudanko/storage/data/hg19/all/ENCODE_tf_peak_calls/wgEncodeRegTfbsClusteredWithCellsV3.bed.gz"
file.K562.DHS <- "/workdir/tc532/manuscript_code/database/merged/min150_sorted_data/K562.bed";

tbo <- read.table(pipe(paste("zcat ", file.TF.chipseq, " | grep K562 | bedtools intersect -a ", file.G1, " -b - -loj")));
tbo <- tbo[,-c(9,10)];

TF.names <- unique(tbo[,8]);
TF.vec0 <- rep(0, NROW(TF.names));
names(TF.vec0) <- TF.names;

tb.dREG <- read.table(file.G1);

mat.TFs <-  do.call("rbind", mclapply( 1:NROW(tb.dREG), function(i){
#mat.TFs <- do.call("rbind", mclapply( 1:200, function(i){
	idx <- which( tbo[,1]== as.character(tb.dREG[i,1]) & tbo[,2]==tb.dREG[i,2] & tbo[,3]==tb.dREG[i,3] );
	TF.vec <- TF.vec0;
	TF.vec[as.character(unique(tbo[idx, 8]))] <- 1;
	return(TF.vec);
  }, mc.cores=30) );


tbh <- read.table(pipe(paste(" bedtools intersect -a ", file.G1, " -b ",  file.K562.DHS, "-wa")));
DHS.status <- unlist(mclapply(1:NROW(tb.dREG), function(i){
	idx <- which( tbh[,1] == as.character(tb.dREG[i,1]) & tbh[,2]==tb.dREG[i,2] & tbh[,3]==tb.dREG[i,3] );
	return( NROW(idx) == 0 );
  }, mc.cores=30) );

colnames(tb.dREG)<-c("chr", "start","end","score");
dreg.mat <- data.frame( tb.dREG, y=DHS.status, mat.TFs);

save(dreg.mat, file="erna.rdata");


train <- sample(1:NROW(dreg.mat))[1:(NROW(dreg.mat)*fold_cv)]
test  <- rep(TRUE, NROW(dreg.mat));
test[train] <- FALSE;
test <- which(test)

## Now the regression.
df <- dreg.mat[,c(4,5)];
sm <- glm(y~score, family=binomial, data=df[train,])	## Strawman model.
df <- dreg.mat[,-c(1:4)];
tf <- glm(y~., family=binomial, data=df[train,])			## Model with ONLY tfs.

df <- dreg.mat[,c(4,5)];
scores_sm <- predict(sm, df[test,])
df <- dreg.mat[,-c(1:4)];
scores_tf <- predict(tf, df[test,])

require(dREG)
roc_sm <- logreg.roc.calc(dreg.mat$y[test], scores_sm)
roc_tf <- logreg.roc.calc(dreg.mat$y[test], scores_tf)

roc.auc(roc_sm)
#0.754
roc.auc(roc_tf)
#0.878

pdf("roc.curve.pdf")
roc.plot(roc_sm, xlim=c(0,1), ylim=c(0,1), col="dark gray")
par(new = TRUE)
roc.plot(roc_tf, xlim=c(0,1), ylim=c(0,1), col="black")
dev.off();


## Now use bootstrap to test the regression coefficients difference from 0.
require(boot);
bb <- boot(data= dreg.mat[,-c(1:4)], R= 1000, statistic= function(a, i) {
   rc <- sample(c(2:NCOL(a)), 1)## Select a random column to leave out.
   vals <- glm(y~., family=binomial, data=a[i,-1*rc])$coefficients
   if(rc==NCOL(a)) {
    ans  <- c(vals, NA)
   } else {
    ans  <- c(vals[1:(rc-1)], NA, vals[(rc):NROW(vals)])
   }
   names(ans) <- c("Intercept", colnames(a[2:NCOL(a)]))
   return(ans)
 },ncpus=30,parallel="multicore")

std.error <- sapply(1:NROW(bb$t0) , function(x) {sd(bb$t[,x], na.rm=TRUE)})
sig <- sapply(1:NROW(bb$t0) , function(x) {!xor(quantile(bb$t[,x], sl, na.rm=TRUE)>0, quantile(bb$t[,x], 1-sl, na.rm=TRUE)>0)}) # 0.025 0.975

source("https://raw.githubusercontent.com/Danko-Lab/dREG/master/dREG_paper_analyses/train_svm/erna_regression/erna_drawbars.R");

pdf("roc.boot.pdf")

drawBars(bb$t0, std.error, names(bb$t0))
drawBars(bb$t0[sig], std.error[sig], names(bb$t0)[sig])
drawBarsVertical(bb$t0[sig], std.error[sig], names(bb$t0)[sig])

dev.off();







#-------------------


dreg_th <- 0.8
dreg_max_neg<- 0.3   ## MAX of the negative samples.

## Read in DNAse-1 peaks
## V7 is signal ... see here: http://genome.ucsc.edu/FAQ/FAQformat.html#format12
#dnase1 <- read.table("/usr/data/GROseq.parser/hg19/k562/dnase/wgEncodeOpenChromDnaseK562PkV2.narrowPeak.gz")

peaks <- read.table(paste(pre, ".overlaps.tsv", sep=""))
tfs <- read.table(paste(pre, ".colNames.txt", sep="")) #read.table("/usr/data/GROseq.parser/hg19/k562/tf_peaks/tf_names_and_files.txt")
colnames(peaks) <- c("chrom", "start", "end", "name", "peakScore", "dreg", as.character(tfs$V3))



pdf(paste(pre,".dnasedreg.pdf", sep=""))

 ## Divide into train and test sets.
 use <- (peaks$dreg>dreg_th | peaks$dreg<dreg_max_neg | is.na(peaks$dreg)) # & (peaks$H3K4me3 == 0) ## Remove peaks just below threshold.

 train <- sample(c(1:NROW(peaks))[use], sum(use)*fold_cv)
 test  <- rep(TRUE, NROW(peaks)); test[train] <- FALSE; test[!use] <- FALSE; test <- which(test)

 ## Create the data structure.
 ## Combine columns that have the same TF (otherwise get some odd combined results).
 ## AND.
 p_combined <- NULL
 for(i in unique(colnames(peaks)[7:NCOL(peaks)])) {
   if(sum(colnames(peaks) == i)>1) {
     p_combined <- cbind(p_combined, as.integer(rowSums(peaks[,colnames(peaks) == i])>0))
   }
   else {
     p_combined <- cbind(p_combined, peaks[,colnames(peaks) == i])
   }
 }
 colnames(p_combined) <- unique(names(peaks)[7:NCOL(peaks)])

 ## Remove certain columns.
 p_combined <- p_combined[,grep("H3K4me3", colnames(p_combined), invert=TRUE)]
 p_combined <- p_combined[,grep("Me3_Me1_ratio", colnames(p_combined), invert=TRUE)]

 data_df <- cbind(data.frame(y= peaks$dreg>dreg_th, peakScore= peaks$peakScore), p_combined) #peaks[,c(12:NCOL(peaks))])
 data_df$y[is.na(data_df$y)] <- FALSE
 data_tfs <- cbind(data.frame(y= peaks$dreg>dreg_th), p_combined) #peaks[,c(12:NCOL(peaks))])
 data_tfs$y[is.na(data_df$y)] <- FALSE

 ## Now the regression.
 sm <- glm(y~peakScore, family=binomial, data=data_df[train,])	## Strawman model.
 tf <- glm(y~., family=binomial, data=data_tfs[train,])			## Model with ONLY tfs.

# ## Now include interaction terms with LASSO variable selection.
# require(glmnet)
# f <- as.formula(y~(.)^2) ## Include all 1st order interaction terms.
# x <- model.matrix(f, data_tfs[train,])
# y <- as.matrix(data_tfs$y[train])

# tf_l <- glmnet(x= x, y= y, family="binomial", alpha=1) ## alpha=1 is the LASSO penalty; 0 is ridge.
# tf_l$beta
# coef(tf_l)
# sort(tf_l$beta[abs(tf_l$beta[,NCOL(tf_l$beta)])>0,NCOL(tf_l$beta)]) ## All positive terms.
# tf_l$beta[grep(":",rownames(tf_l$beta), invert=TRUE),NCOL(tf_l$beta)] ## Positive single terms.
# tf_l$beta[abs(tf_l$beta[,NCOL(tf_l$beta)])>0 & rownames(tf_l$beta)=="EP300",NCOL(tf_l$beta)] ## Specific TF.

 ## Try partial least squares logistic regression.
# require(CMA)
# require(plsgenomics)
# require(glmnet)
# f <- as.formula(y~.) ## Include all 1st order interaction terms.
# x <- model.matrix(f, data_tfs[train,])
# y <- as.factor(data_tfs$y[train])
#
# pls_tf <- pls_lrCMA(X= x, y= y, models=TRUE)
# ftable(pls_tf)
# plot(pls_tf)

 ## Play with PLS-DA, combining PCA with regression.  Similar to pls_lr, above.
# require(pls)
# pls_tf <- mvr(y~., data=data_tfs[train,], validation = "LOO")
# plot(pls_tf, plottype = c("coefficients"))

 ## How well does the model explain the variation?  How accurately can we predict which DNAse-1 sites will have eRNAs?
# scores_fm <- predict(fm, data_df[test,])
 scores_sm <- predict(sm, data_df[test,])
 scores_tf <- predict(tf, data_df[test,])
# scores_tf_l <- predict(tf_l, model.matrix(f, data_tfs[test,]))
# scores_pls <- predict(pls_tf, data_df[test,])

 require(featureDetector)
# roc_fm <- logreg.roc.calc(data_df$y[test], scores_fm)
 roc_sm <- logreg.roc.calc(data_df$y[test], scores_sm)
 roc_tf <- logreg.roc.calc(data_df$y[test], scores_tf)
# roc_tf_l <- logreg.roc.calc(data_df$y[test], scores_tf_l)
# roc_pls<- logreg.roc.calc(data_df$y[test], scores_pls)

# roc.auc(roc_fm)
 roc.auc(roc_sm)
 roc.auc(roc_tf)
# roc.auc(roc_tf_l)
# roc.auc(roc_pls)

 roc.plot(roc_sm, xlim=c(0,1), ylim=c(0,1), col="dark gray")
 par(new = TRUE)
 roc.plot(roc_tf, xlim=c(0,1), ylim=c(0,1), col="black")
# par(new = TRUE)
# roc.plot(roc_fm, xlim=c(0,1), ylim=c(0,1), col="dark red")
# par(new = TRUE)
# roc.plot(roc_tf_l, xlim=c(0,1), ylim=c(0,1), col="green")
# roc.plot(roc_pls, xlim=c(0,1), ylim=c(0,1), col="green")

 ## Now use bootstrap to test the regression coefficients difference from 0.
 require(boot)
 bb <- boot(data= data_tfs[use,], R= 1000, statistic= function(a, i) {
   rc <- sample(c(2:NCOL(a)), 1)## Select a random column to leave out.
   vals <- glm(y~., family=binomial, data=a[i,-1*rc])$coefficients
   if(rc==NCOL(a)) {
    ans  <- c(vals, NA)
   } else {
    ans  <- c(vals[1:(rc-1)], NA, vals[(rc):NROW(vals)])
   }
   names(ans) <- c("Intercept", colnames(a[2:NCOL(a)]))
   return(ans)
 })


save.image(paste(pre,".erna_regression.RData", sep=""))

 std.error <- sapply(1:NROW(bb$t0) , function(x) {sd(bb$t[,x], na.rm=TRUE)})
 sig <- sapply(1:NROW(bb$t0) , function(x) {!xor(quantile(bb$t[,x], sl, na.rm=TRUE)>0, quantile(bb$t[,x], 1-sl, na.rm=TRUE)>0)}) # 0.025 0.975

 source("erna_drawbars.R")
 drawBars(bb$t0, std.error, names(bb$t0))
 drawBars(bb$t0[sig], std.error[sig], names(bb$t0)[sig])
 drawBarsVertical(bb$t0[sig], std.error[sig], names(bb$t0)[sig])

dev.off()

