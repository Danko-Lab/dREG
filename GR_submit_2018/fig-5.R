getProfile <- function( bed1, hMarkFile, dist= 5000, step=25 )
{
	## Load mark.
	hMark <- load.bigWig(hMarkFile)
	mat1 <- NULL;

	if(NROW(bed1)!=0)
	{
		if(NROW(bed1)<100)
		{
			hCountMatrix <- bed.step.bpQuery.bigWig(hMark, center.bed(bed1[,c(1,5,5)], dist, dist), step=step, abs.value=TRUE)
			hmat1 <- matrix(unlist(hCountMatrix), nrow= NROW(bed1), byrow=TRUE);
			mat1 <- colMeans(hmat1);
		}
		else
		{
			hmat1 <- metaprofile.bigWig(center.bed(bed1[,c(1,5,5)], dist, dist), hMark, step=step)
			mat1 <- abs(hmat1$middle)
		}
	}

    unload.bigWig(hMark);

    return(mat1)
}


writeMatplot<- function(bed1, bed2, file.plus.bw, file.minus.bw, hMarkFile, title, subs= NULL, breaks= NULL,
							        cols= NULL, dist= 5000, step=25)
{

	mat1 <- getProfile( bed1, hMarkFile, dist, step )
	mat2 <- getProfile( bed2, hMarkFile, dist, step )

	par( mar=c(5,4,2,2), plt=c(0.15, 0.99, 0.2, 0.8), mgp=c(1.5,0.4,0) );

	if(!is.null(mat1) || !is.null(mat2))
		y.max <- max(c(mat1, mat2))
	else
		y.max <- 10;

	plot(NA, NA, type="n", col="gray", xlim=c(1,ifelse(is.null(mat1), NROW(mat2), NROW(mat1))), ylim=c(0, y.max)/1000,  xlab="Distance (Kbp) ", ylab="", cex.axis=1/3, cex.lab=1/3,  cex.main=1/3, xaxt = "n", yaxt="n", main=title)

	if(!is.null(mat1)) lines(1:NROW(mat1), mat1/1000, col="blue", lwd=1/8);
	if(!is.null(mat2)) lines(1:NROW(mat2), mat2/1000, col="black", lwd=1/8);

	axis(1, c(0, 100, 200, 300, 400), c(-dist/1000, -dist/2000, 0, dist/2000, dist/1000), cex.axis=1/3, cex=1/3 )
	if( y.max>10 )
		axis(2, c(0, y.max/2000, y.max/1000), c(0, round(y.max/2), round(y.max)), cex.axis=1/3, cex=1/3 )
	else if( y.max>=1 )
		axis(2, c(0, y.max/2000, y.max/1000), c(0, round(y.max/2, digits=1), round(y.max,digits=1)), cex.axis=1/3, cex=1/3 )
	else
		axis(2, c(0, y.max/2000, y.max/1000), c(0, signif(y.max/2,digits=2), signif(y.max,digits=2)), cex.axis=1/3, cex=1/3 )

	return()
}


file.plus.bw       <- "../k562/K562_unt.sort.bed.gz_plus.bw"
file.minus.bw      <- "../k562/K562_unt.sort.bed.gz_minus.bw"

file.MNase.bw      <- "/fs/cbsudanko/storage/data/hg19/k562/sydh_mnase/wgEncodeSydhNsomeK562Sig.bigWig"
file.dnase.bw      <- "/fs/cbsudanko/storage/data/hg19/k562/dnase/wgEncodeOpenChromDnaseK562SigV2.bigWig";
file.H3K27me3.bw   <- "/fs/cbsudanko/storage/data/hg19/k562/histones/wgEncodeBroadHistoneK562H3k27me3StdSig.bigWig"
file.H3K4me3.bw    <- "/fs/cbsudanko/storage/data/hg19/k562/histones/wgEncodeBroadHistoneK562H3k4me3StdSig.bigWig"
file.H3K4me1.bw    <- "/fs/cbsudanko/storage/data/hg19/k562/histones/wgEncodeBroadHistoneK562H3k4me1StdSig.bigWig"
file.H3K27ac.bw    <- "/fs/cbsudanko/storage/data/hg19/k562/histones/wgEncodeBroadHistoneK562H3k27acStdSig.bigWig"
file.H3K9me3.bw    <- "/fs/cbsudanko/storage/data/hg19/k562/histones/wgEncodeBroadHistoneK562H3k9me3StdSig.bigWig"

file.dnase.peak    <- "/fs/cbsudanko/storage/data/hg19/k562/dnase/wgEncodeOpenChromDnaseK562PkV2.narrowPeak.gz";
file.H3K27me3.peak <- "/fs/cbsudanko/storage/data/hg19/k562/histones/wgEncodeBroadHistoneK562H3k27me3StdPk.broadPeak.gz"
file.H3K4me3.peak  <- "/fs/cbsudanko/storage/data/hg19/k562/histones/wgEncodeBroadHistoneK562H3k4me3StdPk.broadPeak.gz"
file.H3K4me1.peak  <- "/fs/cbsudanko/storage/data/hg19/k562/histones/wgEncodeBroadHistoneK562H3k4me1StdAln.bed.gz"
file.H3K27ac.peak  <- "/fs/cbsudanko/storage/data/hg19/k562/histones/wgEncodeBroadHistoneK562H3k27acStdAln.bed.gz"

file.G1 <- "../new-rf-201803/G1/G1.dREG.peak.full.bed.gz"

file.TF.chipseq <- "/fs/cbsudanko/storage/data/hg19/all/ENCODE_tf_peak_calls/wgEncodeRegTfbsClusteredWithCellsV3.bed.gz"
file.DHS.Duke <- "/fs/cbsudanko/storage/data/hg19/k562/dnase/wgEncodeOpenChromDnaseK562PkV2.narrowPeak.gz"
file.DHS.UW   <- "/fs/cbsudanko/storage/data/hg19/k562/dnase/GSM646567_hg19_wgEncodeUwDgfK562Pk.narrowPeak.txt.gz"


source("https://raw.githubusercontent.com/Danko-Lab/dREG/master/dREG_paper_analyses/train_svm/erna_regression/erna_drawbars.R");

library(parallel)
library(vioplot)
library(boot)
library(bigWig);
require(dREG);
library(bigWig);


	fold_cv <- 0.9
	sl <- (1e-3)/2

	get_DT<-function()
	{
		tb.dt <- read.table(pipe(paste("zcat ", file.TF.chipseq, " | grep K562 | cut -f 1,2,3,4 | bedtools intersect -a ", file.G1, " -b - -loj")));
		tb.dt <- tb.dt[tb.dt$V7!=".",]
		tb.dt <- unique(tb.dt[,c(1:3,4,6,10)]);

		return(tb.dt);
	}

	tbo <- get_DT();
	tb.dREG <- read.table(file.G1);
	colnames(tb.dREG)<-c("chr", "start","end","score", "prob", "center");

	TF.names <- unique(tbo[,6]);
	TF.vec0 <- rep(0, NROW(TF.names));
	names(TF.vec0) <- TF.names;

	mat.TFs <-  do.call("rbind", mclapply( 1:NROW(tb.dREG), function(i){
	#mat.TFs <- do.call("rbind", mclapply( 1:200, function(i){
		idx <- which( tbo[,1]== as.character(tb.dREG[i,1]) & tbo[,2]==tb.dREG[i,2] & tbo[,3]==tb.dREG[i,3] );
		TF.vec <- TF.vec0;
		TF.vec[as.character(unique(tbo[idx, 6]))] <- 1;
		return(TF.vec);
	  }, mc.cores=30) );


	tbh <- read.table(pipe(paste("bedtools intersect -a ", file.G1, " -b ", file.DHS.Duke,  "-loj")));
	tbh <- tbh[as.character(tbh[,7])!=".", c(1,2,3) ]
	DHS.duke <- unlist(mclapply(1:NROW(tb.dREG), function(i){
		idx <- which( tbh[,1] == as.character(tb.dREG[i,1]) & tbh[,2]==tb.dREG[i,2] & tbh[,3]==tb.dREG[i,3] );
		return( NROW(idx) > 0 );
	  }, mc.cores=30) );


	tbh <- read.table(pipe(paste("bedtools intersect -a ", file.G1, " -b ", file.DHS.UW,  "-loj")));
	tbh <- tbh[as.character(tbh[,7])!=".", c(1,2,3) ]
	DHS.uw <- unlist(mclapply(1:NROW(tb.dREG), function(i){
		idx <- which( tbh[,1] == as.character(tb.dREG[i,1]) & tbh[,2]==tb.dREG[i,2] & tbh[,3]==tb.dREG[i,3] );
		return( NROW(idx) > 0 );
	  }, mc.cores=30) );


	bw.plus <- load.bigWig(file.plus.bw);
	bw.minus <- load.bigWig(file.minus.bw)
	rc.plus <- bed.region.bpQuery.bigWig(bw.plus, tb.dREG[,c(1:3)], op = "sum", abs.value = TRUE);
	rc.minus <- bed.region.bpQuery.bigWig(bw.minus, tb.dREG[,c(1:3)], op = "sum", abs.value = TRUE);
	unload.bigWig(bw.plus);
	unload.bigWig(bw.minus);

	# 2 data models:
	# 0 1 1 1
	#dreg.mat <- data.frame( tb.dREG[,c(1,2,3,6)], sel = TRUE,  y=!( DHS.uw | DHS.duke ), score=tb.dREG[,4], read=rc.plus+rc.minus, mat.TFs);
	# 0 0 0 1
	dreg.mat <- data.frame( tb.dREG[,c(1,2,3,6)], sel = (DHS.duke==DHS.uw),   y=!(DHS.uw & DHS.duke), score=tb.dREG[,4], read=rc.plus+rc.minus, mat.TFs);

	# remove the sites which only intersect with DUKE or UW.
	dreg.mat <- dreg.mat[dreg.mat$sel, -5];

	# remove unbinding peaks:
	dreg.mat <- dreg.mat[-which(rowSums(dreg.mat[,8:NCOL(dreg.mat)])==0), ]

	train <- sample(1:NROW(dreg.mat))[1:(NROW(dreg.mat)*fold_cv)]
	test  <- rep(TRUE, NROW(dreg.mat));
	test[train] <- FALSE;
	test <- which(test)

	## Now the regression.
	df <- dreg.mat[,c(5,6)];
	sm <- glm(y~score, family=binomial, data=df[train,])	## y=score
	df <- dreg.mat[,-c(1:4)];
	tf <- glm(y~., family=binomial, data=df[train,])		## y=score+read+TFs

	df <- dreg.mat[,c(5,6)];
	scores_sm <- predict(sm, df[test,])
	df <- dreg.mat[,-c(1:4)];
	scores_tf <- predict(tf, df[test,])

	roc_sm <- logreg.roc.calc(dreg.mat$y[test], scores_sm);
	roc_tf <- logreg.roc.calc(dreg.mat$y[test], scores_tf);

	roc.auc(roc_sm)
	# 0.8127694
	roc.auc(roc_tf)
	#0.9604431

	pdf("roc.curve.pdf")
	roc.plot(roc_sm, xlim=c(0,1), ylim=c(0,1), col="dark gray")
	par(new = TRUE)
	roc.plot(roc_tf, xlim=c(0,1), ylim=c(0,1), col="black")
	dev.off();

	bb <- boot(data= dreg.mat[,-c(1:4)], R= 1000, statistic= function(a, i) {
	   vals <- glm(y~., family=binomial, data=a[i,])$coefficients
		vals
	},ncpus=30,parallel="multicore")


	tf.ratio <- colSums(dreg.mat[dreg.mat$y,-c(1:7)])/ colSums(dreg.mat[dreg.mat$y,-c(1:7)]);

	pdf("TF.ratio.0001.pdf")
	hist(tf.ratio, breaks=50)
	dev.off();

	std.error <- sapply(1:NROW(bb$t0) , function(x) {sd(bb$t[,x], na.rm=TRUE)})
	sig <- sapply(1:NROW(bb$t0) , function(x) {!xor(quantile(bb$t[,x], sl, na.rm=TRUE)>0, quantile(bb$t[,x], 1-sl, na.rm=TRUE)>0)}) # 0.025 0.975

	save(tf.ratio, tf, dreg.mat, bb, file="roc.0001.boot.rdata");

	pdf("roc.all.0001.boot.pdf")
	drawBars(bb$t0, std.error, names(bb$t0))
	drawBars(bb$t0[sig], std.error[sig], names(bb$t0)[sig])
	drawBarsVertical(bb$t0[sig], std.error[sig], names(bb$t0)[sig])
	dev.off();

}

if(0)
{

	BWs=c(file.plus.bw,
		  file.minus.bw,
		  file.dnase.bw,
		  file.MNase.bw,
		  file.H3K27ac.bw,
		  file.H3K27me3.bw,
		  file.H3K4me1.bw,
		  file.H3K4me3.bw)

	TFs=names(sort(tf.ratio, decreasing=T)[1:20])

	pdf("metaplot.dnase.pdf", width=12, height=6)
	par(mfrow=c(NROW(BWs),NROW(TFs)))
	for(k in 1:NROW(BWs))
	for(i in 1:NROW(TFs))
	{
		tf.status <- dreg.mat[,c(TFs[i])];

		dtdminus <- dreg.mat[dreg.mat$y & tf.status==1, c(1,2,3,4,4)]
		dtd <- dreg.mat[!dreg.mat$y & tf.status==1, c(1,2,3,4,4)]

		title0 <- paste(TFs[i],":",NROW(dtdminus), "/", NROW(dtd), sep="");
		writeMatplot(dtdminus, dtd, file.plus.bw, file.minus.bw, BWs[k],  title=title0, subs= NULL, breaks= NULL, cols= NULL, dist= 5000, step=25)
	}
	dev.off();


	TFs=names(sort(tf.ratio, decreasing=F)[1:20])
	pdf("metaplot.control.pdf", width=12, height=6)
	par(mfrow=c(NROW(BWs),NROW(TFs)))
	for(k in 1:NROW(BWs))
	for(i in 1:NROW(TFs))
	{
		tf.status <- dreg.mat[,c(TFs[i])];

		dtdminus <- dreg.mat[dreg.mat$y & tf.status==1, c(1,2,3,4,4)]
		dtd <- dreg.mat[!dreg.mat$y & tf.status==1, c(1,2,3,4,4)]

		title0 <- paste(TFs[i],":",NROW(dtdminus), "/", NROW(dtd), sep="");
		writeMatplot(dtdminus, dtd, file.plus.bw, file.minus.bw, BWs[k],  title=title0, subs= NULL, breaks= NULL, cols= NULL, dist= 5000, step=25)
	}
	dev.off();


	hist(tf.ratio, breaks=40, xlim=c(0,0.4))
	b <- hist(tf.ratio, breaks=40, plot = F)

	tf.imp <- tf.ratio[c("RPC155", "BRF2", "KAP1", "CHD1", "CEBPB", "NFYB", "GATA2", "TAF7", "GTF3C2", "POLR2A")]
	points( tf.imp, rep(0.5, NROW(tf.imp)), pch=19, cex=0.5, col="blue");

	pos <- floor(tf.imp*100);

	for(i in unique(pos))
	{
		if(NROW(which(pos==i))>1)
			label<-paste( names(tf.imp)[which(pos==i)], collapse=",")
		 else
		 	label<-names(tf.imp)[which(pos==i)];
cat(i,NROW(which(pos==i)), label, "\n");
		text( i/100+0.002, b$counts[i+1]+0.2, label, srt=90, cex=0.75, adj=c(0,1), col="blue" );
	}

	dev.off()


}


if(0)
{
	library(rtfbsdb);

	file.twoBit <- "/fs/cbsudanko/storage/data/hg19/hg19.2bit"

	db <- CisBP.extdata("Homo_sapiens");

	tfs <- tfbs.createFromCisBP(db);

    r.comp <- tfbs.enrichmentTest( tfs, file.twoBit_path, dreg.mat[dreg.mat$y, 1:3], dreg.mat[!dreg.mat$y, 1:3], gc.correction=T, use.cluster=F, ncores = 16, threshold=7, gc.robust.rep=10, background.length=2500);

	tfbs.reportEnrichment(tfs, r.comp, file.pdf = "tfbs.report.pdf", enrichment.type="enriched")

    r.scan <- tfbs.scanTFsite(tfs, file.twoBit, dreg.mat[dreg.mat$y, 1:3], ncores = 7);

	#
    #CEBPB: 3141
    #GATA2: 973,
    #NFYB: 1950,1951
    #SPI1: 757

	TF.idx <- c(313, 973, 1951, 757);
	TF.names <- c("CEBPB","GATA2","NFYB","SPI1");

	BWs=c(file.dnase.bw,
		  file.MNase.bw)

	pdf("rtfbsdb.dnase.pdf")
	par(mfrow=c(NROW(BWs),NROW(TF.idx)))
	for(k in 1:NROW(BWs))
	for(i in 1:4)
	{
		title0 <- TF.names[i];
		bed <- r.scan$result[[ TF.idx[i] ]]
		bed <- bed [bed[,5]>=7,]
		bed[,5] <- round((bed[,2] + bed[,3])/2);
		writeMatplot( bed, NULL, file.plus.bw, file.minus.bw, BWs[k],  title=title0, subs= NULL, breaks= NULL, cols= NULL, dist= 2500, step=25)
	}
	dev.off();




}

if(1)
{

	BWs=c(file.dnase.bw,
		  file.MNase.bw,
		  file.H3K4me1.bw,
		  file.H3K4me3.bw,
		  file.H3K27ac.bw,
		  file.H3K27me3.bw,
		  file.H3K9me3.bw)

	TFs=c("MAZ", "ZNF143", "GATA2", "SPI1", "NFYB", "CEBPB");

	dat.pro <- list();
	for(k in 1:NROW(BWs))
	for(i in 1:NROW(TFs))
	{
		tf.status <- dreg.mat[,c(TFs[i])];
		dtdminus <- dreg.mat[dreg.mat$y & tf.status==1, c(1,2,3,4,4)]
		dtd <- dreg.mat[!dreg.mat$y & tf.status==1, c(1,2,3,4,4)]
		dat.pro[[ k*NROW(TFs) + i ]] <- list( TF=TFs[i], mat1=getProfile( dtdminus, BWs[k], dist= 500, step=5 ), mat2=getProfile( dtd, BWs[k], dist= 500, step=5 ) )
	}

	pdf("metaplot.final2.pdf", width=6, height=6)

#	par( mar=c(5,4,2,2), plt=c(0.15, 0.99, 0.2, 0.8), mgp=c(1.5,0.4,0) );
	par(mfrow = c(NROW(BWs), NROW(TFs)))
	par(cex = 0.6)
	par(mar = c(0, 0, 0, 0), oma = c(4, 4, 1.5, 0.5))
	par(tcl = -0.25)
	par(mgp = c(2, 0.6, 0))

	for(k in 1:NROW(BWs))
	{
		ymax <- 0;
		for(i in 1:NROW(TFs))
		{
			mat1 <- dat.pro[[ k*NROW(TFs) + i ]]$mat1;
			mat2 <- dat.pro[[ k*NROW(TFs) + i ]]$mat2;
			ymax <- max( c(mat1, mat2, ymax ) );
		}

		for(i in 1:NROW(TFs))
		{
			mat1 <- dat.pro[[ k*NROW(TFs) + i ]]$mat1;
			mat2 <- dat.pro[[ k*NROW(TFs) + i ]]$mat2;

			plot(1, axes = FALSE, type = "n", ylim=c(0, ymax), xlim=c(-NROW(mat1)/8, NROW(mat1)*9/8));

			if(k==2)
			{
				segments(-73.5/5 + NROW(mat1)/2, 0, -73.5/5 + NROW(mat1)/2, ymax, col="#888888", lwd=1/4, lty=22);
				segments(+73.5/5 + NROW(mat1)/2, 0, +73.5/5 + NROW(mat1)/2, ymax, col="#888888", lwd=1/4, lty=22);
			}
			if(!is.null(mat1)) lines(1:NROW(mat1), mat1, col="#2ca25f", lwd=0.7);
			if(!is.null(mat2)) lines(1:NROW(mat2), mat2, col="#8856a7", lwd=0.7);

			if (k==NROW(BWs))
				axis(1, col = "grey40", col.axis = "grey20", at =c(0, NROW(mat1)/2, NROW(mat1)), label=c(-5,0,5), cex=0.9)

			if (i==1)
		        axis(2, col = "grey40", col.axis = "grey20", at =c(0, ymax*3/8, ymax*3/4), label=round(c(0, ymax*3/8, ymax*3/4)), cex=0.9)

			box(col = "grey60")

			if(k==1 && i==NROW(TFs))
				legend("topright", legend = c("DNase-", "DNase+"), col=c("#2ca25f", "#8856a7"), lty=1, lwd=0.7, cex=0.9, bty="n");

		}
	}

	for(i in 1:NROW(TFs))
    	mtext(TFs[i], side = 3, outer = TRUE, cex = 0.7, line = 0.2,  at=(i-0.5)/NROW(TFs), col = "grey20")

	BWnames <- c( "H3K9me3", "H3K27me3", "H3K27ac", "H3K4me3", "H3K4me1", "MNase", "DNase");
	for(i in 1:NROW(BWnames))
    	mtext( BWnames[i] , side = 2, outer = TRUE, cex = 0.7, line = 2.2,  at=(i-0.5)/NROW(BWnames), col = "grey20")

   	mtext("Distance(Kb)", side = 1, outer = TRUE, cex = 0.7, line = 2.2,  col = "grey20")

	dev.off();

}


if(1)
{

	BWs=c(file.dnase.bw,
		  file.MNase.bw,
		  file.H3K4me1.bw,
		  file.H3K4me3.bw,
		  file.H3K27ac.bw,
		  file.H3K27me3.bw,
		  file.H3K9me3.bw)

	TFs=c("MAZ", "ZNF143", "GATA2", "SPI1", "NFYB", "CEBPB");

	dat.pro <- list();
	for(k in 1:NROW(BWs))
	for(i in 1:NROW(TFs))
	{
		tf.status <- dreg.mat[,c(TFs[i])];
		dtdminus <- dreg.mat[dreg.mat$y & tf.status==1, c(1,2,3,4,4)]
		dtd <- dreg.mat[!dreg.mat$y & tf.status==1, c(1,2,3,4,4)]
		dat.pro[[ k*NROW(TFs) + i ]] <- list( TF=TFs[i], mat1=getProfile( dtdminus, BWs[k] ), mat2=getProfile( dtd, BWs[k] ) )
	}

	pdf("metaplot.final2.pdf", width=6, height=6)

#	par( mar=c(5,4,2,2), plt=c(0.15, 0.99, 0.2, 0.8), mgp=c(1.5,0.4,0) );
	par(mfrow = c(NROW(BWs), NROW(TFs)))
	par(cex = 0.6)
	par(mar = c(0, 0, 0, 0), oma = c(4, 4, 1.5, 0.5))
	par(tcl = -0.25)
	par(mgp = c(2, 0.6, 0))

	for(k in 1:NROW(BWs))
	{
		ymax <- 0;
		for(i in 1:NROW(TFs))
		{
			mat1 <- dat.pro[[ k*NROW(TFs) + i ]]$mat1;
			mat2 <- dat.pro[[ k*NROW(TFs) + i ]]$mat2;
			ymax <- max( c(mat1, mat2, ymax ) );
		}

		for(i in 1:NROW(TFs))
		{
			mat1 <- dat.pro[[ k*NROW(TFs) + i ]]$mat1;
			mat2 <- dat.pro[[ k*NROW(TFs) + i ]]$mat2;

			plot(1, axes = FALSE, type = "n", ylim=c(0, ymax), xlim=c(-NROW(mat1)/8, NROW(mat1)*9/8));

			if(k==2)
			{
				segments(-73.5/25 + NROW(mat1)/2, 0, -73.5/25 + NROW(mat1)/2, ymax, col="#888888", lwd=1/4, lty=22);
				segments(+73.5/25 + NROW(mat1)/2, 0, +73.5/25 + NROW(mat1)/2, ymax, col="#888888", lwd=1/4, lty=22);
			}
			if(!is.null(mat1)) lines(1:NROW(mat1), mat1, col="#2ca25f", lwd=0.7);
			if(!is.null(mat2)) lines(1:NROW(mat2), mat2, col="#8856a7", lwd=0.7);

			if (k==NROW(BWs))
				axis(1, col = "grey40", col.axis = "grey20", at =c(0, NROW(mat1)/2, NROW(mat1)), label=c(-5,0,5), cex=0.9)

			if (i==1)
		        axis(2, col = "grey40", col.axis = "grey20", at =c(0, ymax*3/8, ymax*3/4), label=round(c(0, ymax*3/8, ymax*3/4)), cex=0.9)

			box(col = "grey60")

			if(k==1 && i==NROW(TFs))
				legend("topright", legend = c("DNase-", "DNase+"), col=c("#2ca25f", "#8856a7"), lty=1, lwd=0.7, cex=0.9, bty="n");

		}
	}

	for(i in 1:NROW(TFs))
    	mtext(TFs[i], side = 3, outer = TRUE, cex = 0.7, line = 0.2,  at=(i-0.5)/NROW(TFs), col = "grey20")

	BWnames <- c( "H3K9me3", "H3K27me3", "H3K27ac", "H3K4me3", "H3K4me1", "MNase", "DNase");
	for(i in 1:NROW(BWnames))
    	mtext( BWnames[i] , side = 2, outer = TRUE, cex = 0.7, line = 2.2,  at=(i-0.5)/NROW(BWnames), col = "grey20")

   	mtext("Distance(Kb)", side = 1, outer = TRUE, cex = 0.7, line = 2.2,  col = "grey20")

	dev.off();

}