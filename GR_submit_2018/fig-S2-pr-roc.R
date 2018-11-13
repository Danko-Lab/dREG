library(dREG)
source("/home/zw355/proj/prj10-dreg/script5/auc-pr.R");

eval_status <- function( positive.bed, tb.negative, pred_infp, sample.ratio=0.1, return.only.label=T)
{
	if( sample.ratio != 1)
	   pred_infp <- pred_infp[ sample(1:NROW(pred_infp))[1:round(NROW(pred_infp)*sample.ratio)], ]

	all_feat <- feat( seqname= pred_infp[,1], start= pred_infp[,2], end= (pred_infp[,3]) )
	positive_feat <- feat(seqname= positive.bed[,1], start= positive.bed[,2], end= positive.bed[,3])
	ol <- overlap.feat(x= all_feat, filter = positive_feat)
	pos_indx <- match(paste(ol$seqname, ol$start, ol$end), paste(all_feat$seqname, all_feat$start, all_feat$end))

	all_feat <- feat( seqname= pred_infp[,1], start= pred_infp[,2], end= (pred_infp[,3]) )
	negative_feat <- feat(seqname= tb.negative[,1], start= tb.negative[,2], end= tb.negative[,3])
	ol <- overlap.feat(x= all_feat, filter = negative_feat)
	neg_indx <- match(paste(ol$seqname, ol$start, ol$end), paste(all_feat$seqname, all_feat$start, all_feat$end))

	y_true <- rep(NA, NROW(pred_infp));
	y_true[pos_indx] <- 1;
	y_true[neg_indx] <- 0;
	pred_bed <- cbind( pred_infp, "TRUE"=y_true);

	return( pred_bed [!is.na(pred_bed$"TRUE"),] )
}

out_holdout_pdf<-function(pred_bed_list, pred_bed_cex, pred_bed_cols, title, file_noext, pt.PRs, ext="PDF")
{
	if(ext=="PDF")
		pdf(paste(file_noext, ".pdf", sep=""))
	else
	{
		png(paste(file_noext, ".pr.png", sep=""), width=600, height=600)
	}

	AUC1 <- AUC2 <- c();
	cols <- c("blue", "purple", "darkgreen", "orange");
	plot(1, 1, xlim=c(0,1), ylim=c(0,1), type="n", xlab="Recall", ylab="Precision", main=title, cex.axis=1.5, cex.lab=1.5);
	for(i in 1:length(pred_bed_list))
	{
		xy<- rocdf( pred_bed_list[[i]][,4], pred_bed_list[[i]][,5], type="pr") ;
		points( xy[,1], xy[,2], cex=pred_bed_cex[i], col=pred_bed_cols[i] );
		AUC1 <- c(AUC1, round( auc_pr ( pred_bed_list[[i]][,5], pred_bed_list[[i]][,4] ), 3 ) );
	}

	#points(rep(1, length(pt.PRs)), pt.PRs, pch=17, col=c("#c45ca2", "#60a862"), cex=1.5);

	lengend_data[1] <- paste(lengend_data[1], "( AUC:", AUC1[1], ")")
	lengend_data[2] <- paste(lengend_data[2], "( AUC:", AUC1[2], ")")
	lengend_data[3] <- paste(lengend_data[3], "( AUC:", AUC1[3], ")")
	lengend_data[4] <- paste(lengend_data[4], "( AUC:", AUC1[4], ")")

	legend("bottomleft",  legend=lengend_data, lwd = lengend_lwd, lty = lengend_lty, pch = lengend_pch, merge = TRUE,
		    col = lengend_cols, text.col = lengend_cols,  pt.cex=lengend_cex, cex=1.5 )

	cat("PR AUC=", AUC1, "\n");

	if(ext!="PDF")
	{
		dev.off();
		png(paste(file_noext, ".roc.png", sep=""), width=800, height=800);
	}


	plot(1, 1, xlim=c(0,1), ylim=c(0,1), type="n", xlab="False Positive Rate", ylab="True Positive Rate", main=title, cex.axis=1.5, cex.lab=1.5);
	for(i in 1:length(pred_bed_list))
	{
		xy<- rocdf( pred_bed_list[[i]][,4], pred_bed_list[[i]][,5], type="roc") ;
		points(xy[,1], xy[,2], cex=pred_bed_cex[i], col=pred_bed_cols[i], pch=21);
		AUC2 <- c(AUC2, round( roc.auc( logreg.roc.calc(pred_bed_list[[i]][,5], pred_bed_list[[i]][,4]) ),3) );
	}

	cat("ROC AUC=", AUC2, "\n");
	for(i in 1:NROW(AUC2) )
		lengend_data[i] <- paste(lengend_data[i], "( AUC:", AUC2[i], ")")

	legend("bottomright", legend=lengend_data, col = pred_bed_cols, text.col = pred_bed_cols, pt.cex=c(pred_bed_cex), cex = 1.5 )

	dev.off();
}


library(rphast)

load("../dreg-train/k562.positive.bed.rdata")
load("../dreg-train/k562.negative.bed.rdata")

#ntb4 <- read.table("../new-model-201709/G4/out.dREG.infp.bed.gz")[,c(1:4)]
#otb4 <- read.table("../old-model/test-G4/out.dREG.pred.gz")[,c(1:4)]
#idx <- match( paste(otb4[,1], otb4[,2], sep=":"), paste(ntb4[,1], ntb4[,2], sep=":"));
#ntb4 <- ntb4[idx[!is.na(idx)],]

#ntb4 <- eval_status( positive_bed, negative_bed, ntb4, sample.ratio=0.1 )
#otb4 <- eval_status( positive_bed, negative_bed, otb4, sample.ratio=0.1 )

ntb7 <- read.table("../new-model-201712/G7/G7.dREG.infp.bed.gz")
#ntb7 <- ntb7[ntb7$V5==1, 1:4];
ntb7 <- ntb7[, 1:4];
otb7 <- read.table("../old-model/test-G7/out.dREG.pred.gz")[,c(1:4)]
idx <- match(paste(ntb7[,1], ntb7[,2], sep=":"), paste(otb7[,1], otb7[,2], sep=":"));
ntb7 <- ntb7[idx[!is.na(idx)],]

ntb7 <- eval_status( positive_bed, negative_bed, ntb7, sample.ratio=1 )
otb7 <- eval_status( positive_bed, negative_bed, otb7, sample.ratio=1 )

ntb8 <- read.table("../new-rf-201803/G8/G8.dREG.infp.bed.gz")
#ntb7 <- ntb7[ntb7$V5==1, 1:4];
ntb8 <- ntb8[, 1:4];
otb8 <- read.table("../old-model/test-G8/out.dREG.pred.gz")[,c(1:4)]
idx <- match(paste(ntb8[,1], ntb8[,2], sep=":"), paste(otb8[,1], otb8[,2], sep=":"));
ntb8 <- ntb8[idx[!is.na(idx)],]

ntb8 <- eval_status( positive_bed, negative_bed, ntb8, sample.ratio=1 )
otb8 <- eval_status( positive_bed, negative_bed, otb8, sample.ratio=1 )

load("../GM12878/GM12878.negative.bed.rdata");
load("../GM12878/GM12878.positive.bed.rdata");

ntb.GM <- read.table("../new-model-201712/GM/GM.dREG.infp.bed.gz")
#ntb.GM <- ntb.GM[ntb.GM$V5==1, 1:4];
ntb.GM <- ntb.GM[, 1:4];
otb.GM <- read.table("../old-model/test-GM12878/out.dREG.pred.gz")[,c(1:4)]
idx <- match(paste(ntb.GM[,1], ntb.GM[,2], sep=":"), paste(otb.GM[,1], otb.GM[,2], sep=":"));
ntb.GM <- ntb.GM[idx[!is.na(idx)],]

ntb.GM <- eval_status( positive_bed, negative_bed, ntb.GM, sample.ratio=1 )
otb.GM <- eval_status( positive_bed, negative_bed, otb.GM, sample.ratio=1 )

load("../HCT116/HCT116.negative.bed.rdata");
load("../HCT116/HCT116.positive.bed.rdata");

ntb.GH <- read.table("../new-rf-201803/GH/GH.dREG.infp.bed.gz")
#ntb.GH <- ntb.GH[ntb.GH$V5==1, 1:4];
ntb.GH <- ntb.GH[, 1:4];
otb.GH <- read.table("../old-model/test-GH/out.dREG.pred.gz")[,c(1:4)]
idx <- match(paste(ntb.GH[,1], ntb.GH[,2], sep=":"), paste(otb.GH[,1], otb.GH[,2], sep=":"));
ntb.GH <- ntb.GH[idx[!is.na(idx)],]

ntb.GH <- eval_status( positive_bed, negative_bed, ntb.GH, sample.ratio=1 )
otb.GH <- eval_status( positive_bed, negative_bed, otb.GH, sample.ratio=1 )

if(0)
{
file.K562.Tfit <- "../Tfit/extend-bigwig/G1-1_bidir_predictions.bed"
file.GM12878.Tfit <-  "../Tfit/extend-bigwig/GM-1_bidir_predictions.bed"

	load("../dreg-train/k562.positive.bed.rdata");
	options(scipen =99) # not to use scientific notation when writing out
	file.pos = tempfile(fileext=".bed");
	write.table(positive_bed, file.pos, quote=F, row.names=F, col.names=F, sep="\t")

	get_score<-function(str){
		s0 <- strsplit(str, ",");
		s1 <- strsplit(s0[[1]][1], "\\|")
		return(as.numeric(s1[[1]][2]));
	}

	tb.TRUE <- read.table(pipe(paste("bedtools intersect -a ", file.K562.Tfit, "-b ",  file.pos, "-wa")))[,c(1:3)]
	#tb.K562 <- data.frame(read.table(file.K562.Tfit)[, c(1:3)], unlist(lapply(as.character(read.table(file.K562.Tfit)[,4]), get_score)), 0);
	tb.K562 <- data.frame(read.table(file.K562.Tfit)[, c(1:3)], 1, 0);
	idx.pos <- match( paste(tb.TRUE[,1], ":", tb.TRUE[,2], ":", tb.TRUE[,3], sep=""), paste(tb.K562[,1], ":", tb.K562[,2], ":", tb.K562[,3], sep="") )
	tb.K562[idx.pos, 5] <- 1;
	pr.K562 <- sum(tb.K562[,5])/NROW(tb.K562)


	tb.TRUE <- read.table(pipe(paste("bedtools intersect -a ", file.GM12878.Tfit, "-b ",  file.pos, "-wa")))[,c(1:3)]
	#tb.GM12878 <- data.frame(read.table(file.GM12878.Tfit)[, c(1:3)], unlist(lapply(as.character(read.table(file.GM12878.Tfit)[,4]), get_score)), 0);
	tb.GM12878 <- data.frame(read.table(file.GM12878.Tfit)[, c(1:3)], 1, 0);
	idx.pos <- match( paste(tb.TRUE[,1], ":", tb.TRUE[,2], ":", tb.TRUE[,3], sep=""), paste(tb.GM12878[,1], ":", tb.GM12878[,2], ":", tb.GM12878[,3], sep="") )
	tb.GM12878[idx.pos, 5] <- 1;
	pr.GM12878 <- sum(tb.GM12878[,5])/NROW(tb.GM12878)
}

pred_bed_list <- list(otb7, ntb7, otb8, ntb8, otb.GM, ntb.GM, otb.GH, ntb.GH);
pred_bed_cols <- c("#c45ca2", "#c45ca2", "#86db63", "#86db63", "#7e96f7", "#7e96f7", "#CCAA80", "#CCAA80");
pred_bed_cex <- c(0.25, 1, 0.25, 1, 0.25, 1, 0.25, 1);

#lengend_data  <- c("K562 by old model", "K562 by new model", "K562 by Tfit", "GM12878 by old model", "GM12878 by new model", "GM12878 by Tfit");
#lengend_lwd <- c(NA, NA, NA, NA, NA, NA);
#lengend_lty <- c(NA, NA, NA, NA, NA, NA);
#lengend_pch <- c(19, 19, 17, 19, 19, 17);
#lengend_cols <- c("#c45ca2", "#c45ca2", "#c45ca2", "#60a862", "#60a862", "#60a862");
#lengend_cex <- c(0.5, 1, 1.5, 0.5, 1, 1.5);

lengend_data  <- c("K562 by old model", "K562 by new model", "K562 by old model", "K562 by new model", "GM12878 by old model", "GM12878 by new model", "HCT116 by old model", "HCT116 by new model");
lengend_lwd <- c(NA, NA, NA, NA, NA, NA, NA, NA);
lengend_lty <- c(NA, NA, NA, NA, NA, NA, NA, NA);
lengend_pch <- c(19, 19, 19, 19, 19, 19, 19, 19);
lengend_cols <- c("#c45ca2", "#c45ca2", "#86db63", "#86db63",  "#7e96f7", "#7e96f7", "#CCAA80", "#CCAA80");
lengend_cex <- c(0.5, 1,  0.5, 1, 0.5, 1, 0.5, 1);

#out_holdout_pdf(pred_bed_list, pred_bed_cex, pred_bed_cols, "", "dREG_holdout_pr_roc_full", c(pr.K562, pr.GM12878),"png");

reduce_points<-function( x, y, cex, col, pch=21, lwd=5)
{
	dfx <- cbind(x=x,y=y);
	dfx <- dfx[is.finite(rowSums(dfx)),];
	
	for(i in 0:199)
	{
		dfx0 <- dfx[x < (i+1)*0.001*5 & x>=i*0.001*5,];
		if(NROW(dfx0)==0)
			next;
			
		dfx1 <- c();
		for(k in 1:2)
		{
			outline.idx <- chull(dfx0);
			if(NROW(outline.idx)>0)
			{
				dfx1 <- rbind(dfx1, dfx0[outline.idx,,drop=FALSE]);
				dfx0 <- dfx0[-outline.idx,]
			}
		}

		outline.idx <- chull(dfx0);
	
		if(i>5)
		polygon( x =  dfx0[outline.idx,1],
		         y =  dfx0[outline.idx,2],
		        border = col, col = col, lwd=lwd  );

		points( x = dfx1[,1],
				y = dfx1[,2],
				cex=cex*2/3, col=col, pch=19);
	}
}

out_holdout_pdf2<-function(pred_bed_list, pred_bed_cex, pred_bed_cols, title, file_noext, pt.PRs, ext="PDF")
{
	if(ext=="PDF")
		pdf(paste(file_noext, ".pdf", sep=""))
	else
	{
		png(paste(file_noext, ".pr.png", sep=""), width=600, height=600)
	}

	AUC1 <- AUC2 <- c();
	cols <- c("blue", "purple", "darkgreen", "orange");
	plot(1, 1, xlim=c(0,1), ylim=c(0,1), type="n", xlab="Recall", ylab="Precision", main=title, cex.axis=1.5, cex.lab=1.5);
	for(i in 1:length(pred_bed_list))
	{
		xy<- rocdf( pred_bed_list[[i]][,4], pred_bed_list[[i]][,5], type="pr") ;
		reduce_points( xy[,1], xy[,2], cex=pred_bed_cex[i], col=pred_bed_cols[i], lwd=ifelse(i %in% c(1,3,5,7), 2, 5) );
		AUC1 <- c(AUC1, round( auc_pr ( pred_bed_list[[i]][,5], pred_bed_list[[i]][,4] ), 3 ) );
	}

	#points(rep(1, length(pt.PRs)), pt.PRs, pch=17, col=c("#c45ca2", "#60a862"), cex=1.5);

    for(i in 1:NROW(lengend_data))
	   lengend_data[i] <- paste(lengend_data[i], "( AUC:", AUC1[i], ")")

	legend("bottomleft",  legend=lengend_data, lwd = lengend_lwd, lty = lengend_lty, pch = lengend_pch, merge = TRUE,
		    col = lengend_cols, text.col = lengend_cols,  pt.cex=lengend_cex, cex=1.0 )

	cat("PR AUC=", AUC1, "\n");

	if(ext!="PDF")
	{
		dev.off();
		png(paste(file_noext, ".roc.png", sep=""), width=800, height=800);
	}

if(0)
{
	plot(1, 1, xlim=c(0,1), ylim=c(0,1), type="n", xlab="False Positive Rate", ylab="True Positive Rate", main=title, cex.axis=1.5, cex.lab=1.5);
	for(i in 1:length(pred_bed_list))
	{
		xy<- rocdf( pred_bed_list[[i]][,4], pred_bed_list[[i]][,5], type="roc") ;
		reduce_points(xy[,1], xy[,2], cex=pred_bed_cex[i], col=pred_bed_cols[i], pch=21);
		AUC2 <- c(AUC2, round( roc.auc( logreg.roc.calc(pred_bed_list[[i]][,5], pred_bed_list[[i]][,4]) ),3) );
	}

	cat("ROC AUC=", AUC2, "\n");
	for(i in 1:NROW(AUC2) )
		lengend_data[i] <- paste(lengend_data[i], "( AUC:", AUC2[i], ")")

	legend("bottomright", legend=lengend_data, col = pred_bed_cols, text.col = pred_bed_cols, pt.cex=c(pred_bed_cex), cex = 1.5 )
}
	dev.off();
}

out_holdout_pdf2(pred_bed_list, pred_bed_cex, pred_bed_cols, "", "dREG_holdout_pr_roc_full", NA ,"PDF");

