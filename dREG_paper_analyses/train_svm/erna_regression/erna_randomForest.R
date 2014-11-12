##
##
## Goal: Regress the presence/ absence of eRNA against 
##       a set of transcription factors.  Identify the 
##       set of TFs that drive eRNA expression, and those
##       that do not.
##

args <- commandArgs(trailingOnly=TRUE)
pre   <- args[1]

dreg_th <- 0.81
dreg_max_neg<- 0.3   ## MAX of the negative samples.
fold_cv <- 0.9
sl <- (1e-3)/2

## Read in DNAse-1 peaks
## V7 is signal ... see here: http://genome.ucsc.edu/FAQ/FAQformat.html#format12
#dnase1 <- read.table("/usr/data/GROseq.parser/hg19/k562/dnase/wgEncodeOpenChromDnaseK562PkV2.narrowPeak.gz")

peaks <- read.table(paste(pre, ".overlaps.tsv", sep=""))
tfs <- read.table(paste(pre, ".colNames.txt", sep="")) #read.table("/usr/data/GROseq.parser/hg19/k562/tf_peaks/tf_names_and_files.txt")
colnames(peaks) <- c("chrom", "start", "end", "name", "col", "strand", "peakscore", "pval", "na", "na", "dreg", as.character(tfs$V3))

indx <- !is.nan(peaks$dreg)
cor.test(peaks$dreg[indx], peaks$peakscore[indx], method="pearson")

pdf(paste(pre,".dnasedreg.randomforest.pdf", sep=""))
# plot(peaks$dreg[indx], peaks$peakscore[indx], xlab="dREG", ylab="DNAse-1 score")
# abline(v=dreg_th, col="red")

 ## Divide into train and test sets.
 trans <- peaks$dreg>dreg_th & !is.na(peaks$dreg)
 notrans <- peaks$dreg<dreg_max_neg | is.na(peaks$dreg)

 use <- (trans | notrans) # & (peaks$H3K4me3 == 0) ## Remove peaks just below threshold.
 matched_indx <- c(sample(which(use & trans), sum(use & trans), replace=FALSE), sample(which(use & notrans), sum(use & trans), replace=FALSE)) ## Create a matched set.
 use <- rep(FALSE, NROW(trans)) ## Make the use variable.
 use[matched_indx] <- TRUE
 
 train <- sample(c(1:NROW(peaks))[use], sum(use)*fold_cv)
 test  <- rep(TRUE, NROW(peaks)); test[train] <- FALSE; test[!use] <- FALSE; test <- which(test)

 ## Create the data structure.
 ## Combine columns that have the same TF (otherwise get some odd combined results).
 ## AND.
 p_combined <- NULL
 for(i in unique(colnames(peaks)[12:NCOL(peaks)])) {
   if(sum(colnames(peaks) == i)>1) {
     p_combined <- cbind(p_combined, as.integer(rowSums(peaks[,colnames(peaks) == i])>0))
   }
   else {
     p_combined <- cbind(p_combined, peaks[,colnames(peaks) == i])
   }
 }
 colnames(p_combined) <- unique(names(peaks)[12:NCOL(peaks)])

 ## Remove certain columns.
 p_combined <- p_combined[,grep("H3K4me3", colnames(p_combined), invert=TRUE)]
 p_combined <- p_combined[,grep("Me3_Me1_ratio", colnames(p_combined), invert=TRUE)]

 data_df <- cbind(data.frame(y= peaks$dreg>dreg_th, peakScore= peaks$peakscore), p_combined) #peaks[,c(12:NCOL(peaks))])
 data_df$y[is.na(data_df$y)] <- FALSE
 data_tfs <- cbind(data.frame(y= data_df$y), p_combined) #peaks[,c(12:NCOL(peaks))])
 data_tfs_rf <- data_tfs
 data_tfs_rf$y <- as.factor(data_tfs_rf$y)
 
 ## Clustering/ PCA
## This is cool, but its not really helping!
# require(cluster)
# plot( agnes(t(data_tfs_rf[,c(2:NCOL(data_tfs_rf))])) )
# 
# pc <- princomp(data_tfs_rf[,-1], cor = TRUE)
# par(pty = "s")
# plot(pc$scores[,1], pc$scores[,2],	xlab = "PC1", ylab = "PC2", type = "n", lwd = 2)
# text(pc$scores[,1], pc$scores[,2], labels = abbreviate(colnames(data_tfs_rf[,-1])), cex = 0.7, lwd = 2)
#
# plot(pc$scores[,2], pc$scores[,3],	xlab = "PC2", ylab = "PC3", type = "n", lwd = 2)
# text(pc$scores[,2], pc$scores[,3], labels = abbreviate(colnames(data_tfs_rf[,-1])), cex = 0.7, lwd = 2)
#
# plot(pc$scores[,3], pc$scores[,4],	xlab = "PC3", ylab = "PC4", type = "n", lwd = 2)
# text(pc$scores[,3], pc$scores[,4], labels = abbreviate(colnames(data_tfs_rf[,-1])), cex = 0.7, lwd = 2)
# 
# plot(pc$scores[,4], pc$scores[,5],	xlab = "PC4", ylab = "PC5", type = "n", lwd = 2)
# text(pc$scores[,4], pc$scores[,5], labels = abbreviate(colnames(data_tfs_rf[,-1])), cex = 0.7, lwd = 2)
 
 
 #######################
 ## Now the Random Forest.
 require(randomForest)
 rf_tfs <- randomForest(y~., data=data_tfs_rf[train,], importance=TRUE)
 plot(rf_tfs)
 cors <- sapply(names(rf_tfs$importance[,4]), function(x) {cor(as.numeric(data_tfs_rf$y), data_tfs_rf[[x]], method="spearman")})
 res_ <- cbind(rf_tfs$importance, rho= cors)
 res_[order(res_[,4]),c(1:2,4:5)]
 
 ## For comparison, still train these.
 sm <- glm(y~peakScore, family=binomial, data=data_df[train,])	## Strawman model.
 fm <- glm(y~., family=binomial, data=data_df[train,])			## Full model with TFs.
 tf <- glm(y~., family=binomial, data=data_tfs[train,])			## Model with ONLY tfs.


 ## How well does the model explain the variation?  How accurately can we predict which DNAse-1 sites will have eRNAs?
 scores_fm <- predict(fm, data_df[test,])
 scores_sm <- predict(sm, data_df[test,])
 scores_tf <- predict(tf, data_df[test,])
 scores_rf <- as.numeric(predict(rf_tfs, data_tfs_rf[test,], type="prob")[,"TRUE"]) ## This is a matrix!  Have to find the correct column!

 require(featureDetector)
 roc_fm <- logreg.roc.calc(data_df$y[test], scores_fm)
 roc_sm <- logreg.roc.calc(data_df$y[test], scores_sm)
 roc_tf <- logreg.roc.calc(data_df$y[test], scores_tf)
 roc_rf <- logreg.roc.calc(data_df$y[test], scores_rf)

 roc.auc(roc_fm)
 roc.auc(roc_sm)
 roc.auc(roc_tf)
 roc.auc(roc_rf)

 roc.plot(roc_sm, xlim=c(0,1), ylim=c(0,1), col="dark gray")
 par(new = TRUE) 
 roc.plot(roc_tf, xlim=c(0,1), ylim=c(0,1), col="black")
 par(new = TRUE) 
 roc.plot(roc_fm, xlim=c(0,1), ylim=c(0,1), col="dark red")
 par(new = TRUE) 
 roc.plot(roc_rf, xlim=c(0,1), ylim=c(0,1), col="dark blue")

 ## Now use bootstrap to test the regression coefficients difference from 0.
 imp <- importance(rf_tfs, scale=FALSE)
 importance <- imp[,2]
 std.error <- rep(1, NROW(importance))
 sig <- abs(importance) > 0

 source("erna_drawbars.R")
 drawBars(importance, std.error, names(imp))
 drawBars(importance[sig], std.error[sig], names(imp)[sig])

dev.off()
