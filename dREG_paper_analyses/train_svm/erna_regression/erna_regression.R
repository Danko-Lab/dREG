##
##
## Goal: Regress the presence/ absence of eRNA against 
##       a set of transcription factors.  Identify the 
##       set of TFs that drive eRNA expression, and those
##       that do not.
##

args <- commandArgs(trailingOnly=TRUE)
pre   <- args[1]

dreg_th <- 0.8
dreg_max_neg<- 0.3   ## MAX of the negative samples.
fold_cv <- 0.9
sl <- (1e-3)/2

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

