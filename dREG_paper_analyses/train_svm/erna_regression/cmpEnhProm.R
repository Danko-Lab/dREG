##
##
## Goal: Compare regressions between promoters and enhancers 
##       determined using histone marks.  Asking: Are there 
##       any TFs that differ between promoters and enhancers?
##

args <- commandArgs(trailingOnly=TRUE)
pre   <- args[1]

dreg_th <- 0.7
dreg_max_neg<- 0.4   ## MAX of the negative samples.
fold_cv <- 0.9
sl <- (1e-3)/2

## Read in DNAse-1 peaks
## V7 is signal ... see here: http://genome.ucsc.edu/FAQ/FAQformat.html#format12
peaks <- read.table(paste(pre, ".overlaps.tsv", sep=""))
tfs <- read.table(paste(pre, ".colNames.txt", sep="")) #read.table("/usr/data/GROseq.parser/hg19/k562/tf_peaks/tf_names_and_files.txt")
colnames(peaks) <- c("chrom", "start", "end", "name", "col", "strand", "peakscore", "pval", "na", "na", "dreg", as.character(tfs$V3))

pdf(paste(pre,".dnasedreg.pdf", sep=""))
 hist(peaks$Me3_Me1_ratio, 150)
 hist(peaks$Me3_Me1_ratio[peaks$dreg > dreg_th], 150)
 
 ## Divide into train and test sets.
 enh <- (peaks$Me3_Me1_ratio < 0) & (peaks$dreg>dreg_th | peaks$dreg<dreg_max_neg | is.na(peaks$dreg))
 pro <- (peaks$Me3_Me1_ratio > 0) & (peaks$dreg>dreg_th | peaks$dreg<dreg_max_neg | is.na(peaks$dreg))
 to_rand<- which((peaks$dreg>dreg_th | peaks$dreg<dreg_max_neg | is.na(peaks$dreg)))
 
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

 ## Now the regression. Compare enhancers and promoters.
 tf_e <- glm(y~., family=binomial, data=data_tfs[enh,])
 tf_p <- glm(y~., family=binomial, data=data_tfs[pro,])
 ec <- tf_e$coefficients
 pc <- tf_p$coefficients
 
 ## Get the correlation and compare it to a random grouping.
 cor.test(ec, pc, method="spearman") ## TODO: Is this any different from a random grouping?!

 cor_dist <- NULL
 for(i in 1:10) {
  enh_permut <- rep(FALSE, NROW(to_rand))
  enh_permut[sample(1:NROW(to_rand), sum(enh), replace=FALSE)] <- TRUE
  pro_permut <- rep(TRUE, NROW(to_rand))
  pro_permut[enh_permut] <- FALSE
  
  ec_permut <- glm(y~., family=binomial, data=data_tfs[enh_permut,])$coefficients
  pc_permut <- glm(y~., family=binomial, data=data_tfs[pro_permut,])$coefficients
  cor_dist <- c(cor_dist, cor.test(ec_permut, pc_permut, method="spearman")$estimate[[1]])
 }
 
 ## Plot the data.
 plot(ec, pc, xlim=c(-2,2), ylim=c(-2,2))
 abline(0,1)
 boxplot(peaks[enh,"dreg"], peaks[pro,"dreg"], names=c("Enhancers", "Promoters"))

 ## Now use bootstrap to test the regression coefficients difference from 0.
 require(boot)
 bb_e <- boot(data= data_tfs[enh,], R= 100, statistic= function(a, i) {
   glm(y~., family=binomial, data=a[i,])$coefficients
 })
 bb_p <- boot(data= data_tfs[pro,], R= 100, statistic= function(a, i) {
   glm(y~., family=binomial, data=a[i,])$coefficients
 })

 save.image(paste(pre,".e_p_boot.RData", sep=""))

 # std.error <- sapply(1:NROW(bb$t0) , function(x) {sd(bb$t[,x], na.rm=TRUE)})
 # sig <- sapply(1:NROW(bb$t0) , function(x) {!xor(quantile(bb$t[,x], sl, na.rm=TRUE)>0, quantile(bb$t[,x], 1-sl, na.rm=TRUE)>0)}) # 0.025 0.975

 # source("erna_drawbars.R")
 # drawBars(bb$t0, std.error, names(bb$t0))
 # drawBars(bb$t0[sig], std.error[sig], names(bb$t0)[sig])

dev.off()
