##
## Creates heatmaps on human dREG-HD sites.

require(bigWig)
require(pheatmap)
require(RColorBrewer)
library(gtable)
library(grid)
library(sqldf);

options(scipen =99)

rowMax <- function(x) { sapply(1:NROW(x), function(i) {return(max(x[i,], na.rm=TRUE))}) }

transform.bed <- function(anchor, bed, upstreamWindow, downstreamWindow) {
  start = anchor - upstreamWindow
  end = anchor + downstreamWindow + 1

  N = dim(bed)[2]
  if (N >= 6) { # use strands
    is_minus = bed[,6] == '-'
    start[is_minus] = anchor[is_minus] - downstreamWindow
    end[is_minus] = anchor[is_minus] + upstreamWindow + 1
  }
  start = as.integer(start)
  end = as.integer(end)

  res = NULL
  if (N > 3)
    res = data.frame(bed[,1], start, end, bed[,4:N])
  else
    res = data.frame(bed[,1], start, end)

  colnames(res) <- colnames(bed)
  rownames(res) <- rownames(bed)

  return(res)
}

center.bed <- function(bed, upstreamWindow, downstreamWindow) {
  anchor = bed[,2] + ((bed[,3] - bed[,2]) %/% 2)

  transform.bed(anchor, bed, upstreamWindow, downstreamWindow)
}


dist_X_center <- function( pred.bed, file.X )
{
	options("scipen"=100, "digits"=4);

	tmp.bed <- tempfile(fileext=".bed" )
	write.table( pred.bed[,c(1:3)], file=tmp.bed, row.names=F, col.names=F, quote=F, sep="\t");
	system( paste("sort-bed ",  tmp.bed, " > ", tmp.bed, ".sorted", sep="") );

	tb.X <- read.table( file.X, header=F);
	tb.X[,2] <- round((tb.X[,2] + tb.X[,3])/2)
	tb.X[,3] <- tb.X[,2] + 1;
	tmp2.bed <- tempfile(fileext=".bed")
	write.table(tb.X[, c(1:3)], file=tmp2.bed, row.names=F, col.names=F, quote=F, sep="\t");

	tb.close <-  read.table( file = pipe(paste("sort-bed ",  tmp2.bed, " | bedtools closest -d -a ", tmp.bed, ".sorted -b - -t first", sep="") ) );
	return(tb.close);
}

dist_X_close <- function( pred.bed, file.X )
{
	options("scipen"=100, "digits"=4);

	tmp.bed <- tempfile(fileext=".bed" )
	write.table( pred.bed[,c(1:3)], file=tmp.bed, row.names=F, col.names=F, quote=F, sep="\t");
	system( paste("sort-bed ",  tmp.bed, " > ", tmp.bed, ".sorted", sep="") );

	tb.X <- read.table( file.X, header=F);
	tb.X0 <- tb.X[tb.X$V5=="-",];
	tb.X1 <- tb.X[tb.X$V5=="+",];
	tb.X0$V2 <- tb.X0$V3 - 1;
	tb.X1$V3 <- tb.X1$V2 + 1;

	tmp2.bed <- tempfile(fileext=".bed")
	write.table( rbind(tb.X0[, c(1:3)], tb.X1[, c(1:3)]),  file=tmp2.bed, row.names=F, col.names=F, quote=F, sep="\t");
	tb.close <-  read.table( file = pipe(paste("sort-bed ",  tmp2.bed, " | bedtools closest -d -a ", tmp.bed, ".sorted -b - -t first", sep="") ) );

	return(tb.close);
}


pheatmap_2layer <- function( hmat.plus, hmat.minus, breaks.plus=NULL, breaks.minus=NULL, cluster_rows = FALSE, cluster_cols = FALSE, legend=FALSE, show_rownames=FALSE, show_colnames=FALSE, silent=T )
{
	get_color_map <- function(pMat, mMat)
	{
		cMap <- matrix("", nrow=NROW(pMat), ncol=NCOL(pMat));

		if(is.null(breaks.minus)) breaks.minus <- hmat.minus;
		mMat <-  ( mMat -  min(breaks.minus)) / (max(breaks.minus) -  min(breaks.minus));
		if ( sum(mMat<0)>0 ) mMat[mMat<0] <- 0;
		if ( sum(mMat>1)>0 ) mMat[mMat>1] <- 1;

		if(is.null(breaks.plus)) breaks.plus <- hmat.plus;
		pMat <-  ( pMat -  min(breaks.plus)) / (max(breaks.plus) -  min(breaks.plus));
		if ( sum(pMat<0)>0 ) pMat[pMat<0] <- 0;
		if ( sum(pMat>1)>0 ) pMat[pMat>1] <- 1;

		for(i in 1:NROW(pMat))
		for(j in 1:NCOL(pMat))
		{
			B <- max( mMat[i,j], 1-pMat[i,j] );
			G <- 1 - max( mMat[i,j], pMat[i,j] );
			R <- max( 1- mMat[i,j], pMat[i,j] );

			cMap[i,j] <- rgb(R, G,B);
		}
		return(cMap);
	}

	gt <- pheatmap( hmat.plus, cluster_rows = cluster_rows, cluster_cols = cluster_cols, col= rgb( 0, 0, 1 - breaks.plus/max(breaks.plus)), breaks = breaks.plus, legend=legend, show_rownames=show_rownames, show_colnames=show_colnames, silent=silent );
	gt$gtable$grobs[[1]]$children[[1]]$gp$fill <- get_color_map( hmat.plus, hmat.minus);

	return(gt);
}



writeProseqHeatmap<- function(bed, file.plus.bw, file.minus.bw, png.name, subs= NULL, breaks= NULL,
							cols= NULL, dist= 5000, step=25)
{
	navg <- 20 ## Average every navg rows

	# distance
	tb.close <- dist_X_close(bed[,c(1:3)], file.protein.coding );
	idx.proxm <- which( tb.close$V7 <= 0 )
	idx.distal <- which( tb.close$V7 > 1000)
	cat("<=1K", NROW(idx.proxm), ">1K", NROW(idx.distal), "\n");

	## Load GROseq data .
	hPlus <- load.bigWig(file.plus.bw)
	hMinus <- load.bigWig(file.minus.bw)

	if(0)
	{
	## Get a matrix of counts.
	hPlusMatrix <- bed.step.bpQuery.bigWig(hPlus, center.bed(bed[,c(1,5,5)], dist, dist), step=step, abs.value=TRUE)
	hmat <- log(matrix(unlist(hPlusMatrix), nrow= NROW(bed), byrow=TRUE)+1)
	hm_order   <- order(bed$score, decreasing=TRUE);
	if(is.null(hm_order)) {
	  hm_order <- order(rowSums(hmat[,(NCOL(hmat)/2 -10):(NCOL(hmat)/2 +10)]), decreasing=TRUE)
	}
	hmat.plus <- hmat[hm_order,]

	hMinusMatrix <- bed.step.bpQuery.bigWig(hMinus, center.bed(bed[,c(1,5,5)], dist, dist), step=step, abs.value=TRUE)
	hmat <- log(matrix(unlist(hMinusMatrix), nrow= NROW(bed), byrow=TRUE)+1)
	hm_order   <- order(bed$score, decreasing=TRUE);
	if(is.null(hm_order)) {
	  hm_order <- order(rowSums(hmat[,(NCOL(hmat)/2 -10):(NCOL(hmat)/2 +10)]), decreasing=TRUE)
	}
	hmat.minus <- hmat[hm_order,]

	## Average by rows of 10.
	navg <- 20 ## Average every navg rows
	avgMat <- t(sapply(1:floor(NROW(hmat.plus)/navg), function(x) {colMeans(hmat.plus[((x-1)*navg+1):min(NROW(hmat.plus),(x*navg)),])}))
	hmat.plus <- avgMat

	navg <- 20 ## Average every navg rows
	avgMat <- t(sapply(1:floor(NROW(hmat.minus)/navg), function(x) {colMeans(hmat.minus[((x-1)*navg+1):min(NROW(hmat.minus),(x*navg)),])}))
	hmat.minus <- avgMat
	}

	##proximal
	hCountMatrix <- bed.step.bpQuery.bigWig(hPlus, center.bed(bed[idx.proxm,c(1,5,5)], dist, dist), step=step, abs.value=TRUE)
	hmat <- log(matrix(unlist(hCountMatrix), nrow= NROW(idx.proxm), byrow=TRUE)+1)
	hmat.ord <- hmat[order(bed$score[idx.proxm], decreasing=TRUE),]
	hmat.p.plus <- t(sapply(1:floor(NROW(hmat.ord)/navg), function(x) {colMeans(hmat.ord[((x-1)*navg+1):min(NROW(hmat.ord),(x*navg)),])}))

	hCountMatrix <- bed.step.bpQuery.bigWig(hMinus, center.bed(bed[idx.proxm,c(1,5,5)], dist, dist), step=step, abs.value=TRUE)
	hmat <- log(matrix(unlist(hCountMatrix), nrow= NROW(idx.proxm), byrow=TRUE)+1)
	hmat.ord <- hmat[order(bed$score[idx.proxm], decreasing=TRUE),]
	hmat.p.minus <- t(sapply(1:floor(NROW(hmat.ord)/navg), function(x) {colMeans(hmat.ord[((x-1)*navg+1):min(NROW(hmat.ord),(x*navg)),])}))

	##distal
	hCountMatrix <- bed.step.bpQuery.bigWig(hPlus, center.bed(bed[idx.distal,c(1,5,5)], dist, dist), step=step, abs.value=TRUE)
	hmat <- log(matrix(unlist(hCountMatrix), nrow= NROW(idx.distal), byrow=TRUE)+1)
	hmat.ord <- hmat[order(bed$score[idx.distal], decreasing=TRUE),]
	hmat.d.plus <- t(sapply(1:floor(NROW(hmat.ord)/navg), function(x) {colMeans(hmat.ord[((x-1)*navg+1):min(NROW(hmat.ord),(x*navg)),])}))

	hCountMatrix <- bed.step.bpQuery.bigWig(hMinus, center.bed(bed[idx.distal,c(1,5,5)], dist, dist), step=step, abs.value=TRUE)
	hmat <- log(matrix(unlist(hCountMatrix), nrow= NROW(idx.distal), byrow=TRUE)+1)
	hmat.ord <- hmat[order(bed$score[idx.distal], decreasing=TRUE),]
	hmat.d.minus <- t(sapply(1:floor(NROW(hmat.ord)/navg), function(x) {colMeans(hmat.ord[((x-1)*navg+1):min(NROW(hmat.ord),(x*navg)),])}))

	## Write out a heatmap.
	if(is.null(breaks)) {
		bk <- seq(min(min(hmat.p.plus), min(hmat.d.plus)), max(max(hmat.p.plus), max(hmat.d.plus))+0.01, 0.01)
	} else {
		bk <- breaks
	}

	if(is.null(cols)) {
		hmcols.plus <- rev(colorRampPalette(brewer.pal(9,"RdBu"))(length(bk)-1))
	} else {
		hmcols.plus <- colorRampPalette(cols)(length(bk)-1) # red
	}
	breaks.plus <- bk;

	if(is.null(breaks)) {
		bk <- seq(min(min(hmat.p.minus), min(hmat.d.minus)), max(max(hmat.p.minus), max(hmat.d.minus))+0.01, 0.01)
	} else {
		bk <- breaks
	}

	if(is.null(cols)) {
		hmcols.minus <- rev(colorRampPalette(brewer.pal(9,"RdBu"))(length(bk)-1))
	} else {
		hmcols.minus <- colorRampPalette(cols)(length(bk)-1) # red
	}
	breaks.minus <- bk;


	png(paste(png.name,".png",sep=""), width=450, height = 800*1.55 )

	n.all <- NROW(hmat.d.plus) + NROW(hmat.p.plus);
	lay.heights <- c(NROW(hmat.p.plus)/n.all, NROW(hmat.d.plus)/n.all, 0.4, 0.15);
	lay.widths  <- c(1, 4, 5)
	layout(matrix(c(1,2,2, 3,4,4, 5,5,5, 6,6,7 ), nrow=4, byrow=T), widths=lay.widths, heights=lay.heights)

	##part1:
	par(mar=c(0,0,0,0), plt=c(0.2, 0.8,0.01, 0.99 ));
	plot(NA,NA, type="n", xlab="", ylab="", xlim=c(0,1), ylim=c(0, NROW(hmat.p.plus)), xaxs="i", yaxs="i", xaxt = "n", yaxt = "n", bty="n" );
	text(0.5, NROW(hmat.p.plus)/2, paste("TSS(", NROW(idx.proxm), ")", sep=""), adj=c(0.5, 0.5), srt=90, cex=4.5);

	##part2:
	gt <- pheatmap_2layer( hmat.p.plus, hmat.p.minus, breaks.plus, breaks.minus );

	par(mar=c(0,0,0,0), plt=c(0.2, 0.8,0.2, 0.8 ));
    plot(NA,NA, type="n", xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i", xaxt = "n", yaxt = "n");
    ##grid.newpage()
	pushViewport(viewport(layout = grid.layout(4, 3, widths=lay.widths, heights=lay.heights) ))
	gt$gtable$vp <- viewport(layout.pos.row = 1, layout.pos.col = c(2,3))
    grid.draw(gt$gtable)
	popViewport()

	##part3:
	par(mar=c(0,0,0,0), plt=c(0.2, 0.8,0.01, 0.99 ));
	plot(NA,NA, type="n", xlab="", ylab="", xlim=c(0,1), ylim=c(0, NROW(hmat.d.plus)), xaxs="i", yaxs="i", xaxt = "n", yaxt = "n", bty="n" );
	text(0.5, NROW(hmat.d.plus)/2, paste("Distal(", NROW(idx.distal), ")", sep=""), adj=c(0.5, 0.5), srt=90, cex=4.5);

	##part4:
	gt <- pheatmap_2layer( hmat.d.plus, hmat.d.minus, breaks.plus, breaks.minus );

	par(mar=c(0,0,0,0), plt=c(0.2, 0.8,0.2, 0.8 ));
    plot(NA,NA, type="n", xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i", xaxt = "n", yaxt = "n");
    ##grid.newpage()
	pushViewport(viewport(layout = grid.layout(4, 3, widths=lay.widths, heights=lay.heights)))
	gt$gtable$vp <- viewport(layout.pos.row = 2, layout.pos.col =  c(2,3))
	grid.draw(gt$gtable)
	popViewport()

	##part 5:
	## draw mat-plot for overlapping peak and non-overlapped peak

	mat.p <- colMeans( rbind(2^hmat.d.plus, 2^hmat.p.plus) );
	hCountMatrix <- bed.step.bpQuery.bigWig(hPlus, center.bed(bed[,c(1,5,5)], dist, dist), step=step, abs.value=TRUE)
	mat.p <- colMeans( matrix(unlist(hCountMatrix), nrow= NROW(bed), byrow=TRUE) );

	mat.m <- -1*colMeans( rbind(2^hmat.d.minus, 2^hmat.p.minus) );
	hCountMatrix <- bed.step.bpQuery.bigWig(hMinus, center.bed(bed[,c(1,5,5)], dist, dist), step=step, abs.value=TRUE)
	mat.m <- -1*colMeans( matrix(unlist(hCountMatrix), nrow= NROW(bed), byrow=TRUE) );

	par(mar=c(10,0,0,0), plt=c(0.125, 0.99, 0.2, 0.90), mgp=c(5,2,0) );

	y.max <- max(c(mat.p, -1*mat.m));
	plot(NA, NA, type="n", col="gray", xlim=c(1,NROW(mat.p)), ylim=c(-y.max, y.max)/1000,  xlab="Distance (Kbp) ", ylab="", cex.axis=3, cex.lab=3,  xaxt = "n", yaxt="n")
	lines(1:NROW(mat.m), mat.m/1000, col="red", lwd=2);
	lines(1:NROW(mat.p), mat.p/1000, col="blue", lwd=2);
	axis(1, c(0, 100, 200, 300, 400), c(-dist/1000, -dist/2000, 0, dist/2000, dist/1000), cex.axis=3, cex=3 );
	if( y.max>0 )
		axis(2, c(-y.max/1000, 0, y.max/1000), c(round(-y.max), 0, round(y.max)), cex.axis=3, cex=3 )
	else
		axis(2, c(-y.max/1000, 0, y.max/1000), c(signif(-1*y.max,digits=2), 0, signif(y.max,digits=2)), cex.axis=3, cex=3 )


	##part 6:
	## draw colorScale for Heatmap
	bk <- breaks.minus;
	hmcols <- rgb(1, 1-breaks.minus/max(breaks.minus), 1-breaks.minus/max(breaks.minus) );
	{
		par(mar=c(10,0,0,0), plt=c(0.15, 0.85, 0.4, 0.8 ));
		plot(NA,NA, type="n", xlim=c(1, NROW(bk)), ylim=c(0,1), xlab="", ylab="", bty="n", xaxt="n", yaxt= "n" );
		for(i in 1:NROW(bk)) rect(i,0, i+1, 1, col=hmcols[i], border=hmcols[i]);
		axis(1, c(1, NROW(bk)/2, NROW(bk)), round(exp(c(bk[1], bk[round(NROW(bk)/2)], bk[NROW(bk)])),1), cex.axis=3, cex=3, tick=FALSE );
	}

	##part 7:
	## draw colorScale for peak track
	bk <- breaks.plus;
	hmcols <- rgb(1-breaks.plus/max(breaks.plus), 1-breaks.plus/max(breaks.plus), 1 );
	{
		par(mar=c(10,0,0,0), plt=c(0.15, 0.85, 0.4, 0.8 ));
		plot(NA,NA, type="n", xlim=c(1, NROW(bk)), ylim=c(0,1), xlab="", ylab="", bty="n", xaxt="n", yaxt= "n" );
		for(i in 1:NROW(bk)) rect(i,0, i+1, 1, col=hmcols[i], border=hmcols[i]);
		axis(1, c(1, NROW(bk)/2, NROW(bk)), round(exp(c(bk[1], bk[round(NROW(bk)/2)], bk[NROW(bk)])),1), cex.axis=3, cex=3, tick=FALSE );
	}


	dev.off();

	return()
}

writeColorScale<- function(breaks, cols, name) {
     ## Write out a heatmap.
     if(is.null(breaks)) {
                bk <- seq(min(hmat), max(hmat), 0.01)
      } else {
              bk <- breaks
      }

     if(is.null(cols)) {
             hmcols <- rev(colorRampPalette(brewer.pal(9,"RdBu"))(length(bk)-1))
     } else {
            hmcols <- colorRampPalette(cols)(length(bk)-1) # red
     }

	png(paste(name,".scale.png",sep=""), width=500, height=300)
	pheatmap(rev(bk), cluster_rows = FALSE, cluster_cols = FALSE, col= hmcols, breaks = bk, legend=TRUE, legend_breaks= quantile(bk), legend_labels= signif(exp(quantile(bk)),3), show_rownames=FALSE, show_colnames=FALSE)
	dev.off()
}


drawDregHeatmap<-function(dregX, file.prefix)
{
	dregX <- dregX[dregX$V5<=0.05,];
	dregX <- dregX[,-5];
	colnames(dregX) <- c("chr", "start", "end", "score", "peak")
	dregX <- dregX[grep("_|chrY", dregX[,1], invert=TRUE),]
	dregX <- cbind(dregX, width=dregX$end-dregX$start);

	sup <- writeProseqHeatmap( dregX, file.plus.bw, file.minus.bw, paste(file.prefix, "-Proseq",sep="") )
}

if(1)
{

	#MNase.path="/fs/cbsudanko/storage/data/hg19/k562/sydh_mnase/"
	#file.MNase.bw <- "wgEncodeSydhNsomeK562Sig.bigWig"

	file.grocap        <- ""

	#file.H3K9me3.bw    <- "/fs/cbsudanko/storage/data/hg19/k562/histones/wgEncodeBroadHistoneK562H3k9me3StdSig.bigWig"
	file.H3K27me3.bw   <- "/fs/cbsudanko/storage/data/hg19/k562/histones/wgEncodeBroadHistoneK562H3k27me3StdSig.bigWig"
	file.H3K4me3.bw    <- "/fs/cbsudanko/storage/data/hg19/k562/histones/wgEncodeBroadHistoneK562H3k4me3StdSig.bigWig"
	file.H3K4me1.bw    <- "/fs/cbsudanko/storage/data/hg19/k562/histones/wgEncodeBroadHistoneK562H3k4me1StdSig.bigWig"
	file.H3K27ac.bw    <- "/fs/cbsudanko/storage/data/hg19/k562/histones/wgEncodeBroadHistoneK562H3k27acStdSig.bigWig"

	#file.H3K9me3.peak  <- "/fs/cbsudanko/storage/data/hg19/k562/histones/wgEncodeBroadHistoneK562H3k9me3StdAln.bed.gz"
	file.H3K27me3.peak <- "/fs/cbsudanko/storage/data/hg19/k562/histones/wgEncodeBroadHistoneK562H3k27me3StdPk.broadPeak.gz"
	file.H3K4me3.peak  <- "/fs/cbsudanko/storage/data/hg19/k562/histones/wgEncodeBroadHistoneK562H3k4me3StdPk.broadPeak.gz"
	file.H3K4me1.peak  <- "/fs/cbsudanko/storage/data/hg19/k562/histones/wgEncodeBroadHistoneK562H3k4me1StdAln.bed.gz"
	file.H3K27ac.peak  <- "/fs/cbsudanko/storage/data/hg19/k562/histones/wgEncodeBroadHistoneK562H3k27acStdAln.bed.gz"

	file.dnase.bw = "/fs/cbsudanko/storage/data/hg19/k562/dnase/wgEncodeOpenChromDnaseK562SigV2.bigWig";
	file.dnase.peak = "/fs/cbsudanko/storage/data/hg19/k562/dnase/wgEncodeOpenChromDnaseK562PkV2.narrowPeak.gz";

	file.grocap <- "../k562/hg19.k562.grocap.pair.bed";
	file.protein.coding <- "./gencode.v19.protein.coding.bed"

	#file.plus.bw <- "../k562/6045_7157_27176_HNHKJBGXX_K562_0min_celastrol10uM_rep2_GB_CAGATC_R1_plus.primary.bw"
	#file.minus.bw <- "../k562/6045_7157_27176_HNHKJBGXX_K562_0min_celastrol10uM_rep2_GB_CAGATC_R1_minus.primary.bw"
	file.plus.bw <- "../k562/K562_unt.sort.bed.gz_plus.bw"
	file.minus.bw <- "../k562/K562_unt.sort.bed.gz_minus.bw"

	dregP <- read.table("../new-rf-201803/G1/G1.dREG.peak.full.bed.gz");
	drawDregHeatmap(dregP, "dregG1");
}



