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


writeHeatmap<- function(bed, file.plus.bw, file.minus.bw, hMarkFile, hPeakFile, png.name, subs= NULL, breaks= NULL,
							cols= NULL, dist= 5000, step=25)
{
	## Load mark.
	hMark <- load.bigWig(hMarkFile)

	## Get a matrix of counts.
	hCountMatrix <- bed.step.bpQuery.bigWig(hMark, center.bed(bed[,c(1,5,5)], dist, dist), step=step, abs.value=TRUE)
	hmat <- log(matrix(unlist(hCountMatrix), nrow= NROW(bed), byrow=TRUE)+1)

	hm_order   <- order(bed$score, decreasing=TRUE);
	if(is.null(hm_order)) {
	  hm_order <- order(rowSums(hmat[,(NCOL(hmat)/2 -10):(NCOL(hmat)/2 +10)]), decreasing=TRUE)
	}
	hmat.org <- hmat[hm_order,]

	## Average by rows of 20.
	navg <- 20 ## Average every navg rows
	avgMat <- t(sapply(1:floor(NROW(hmat.org)/navg), function(x) {colMeans(hmat.org[((x-1)*navg+1):min(NROW(hmat.org),(x*navg)),])}))
	hmat.all <- avgMat

	#hCountSum <- bed.region.bpQuery.bigWig( hMark, bed[,c(1,2,3)])
	#hmat.ord <- hCountSum[order(bed$score, decreasing=TRUE)]
	#avgMat.all <- sapply(1:floor(NROW(hmat.ord)/navg), function(x) {mean(hmat.ord[((x-1)*navg+1):min(NROW(hmat.ord),(x*navg))])})

	## Write out a heatmap.
	if(is.null(breaks)) {
		bk <- seq(min(hmat.all), max(hmat.all), 0.01)
	} else {
		bk <- breaks
	}

	if(is.null(cols)) {
		hmcols <- rev(colorRampPalette(brewer.pal(9,"RdBu"))(length(bk)-1))
	} else {
		hmcols <- colorRampPalette(cols)(length(bk)-1) # red
	}

	# distance
	tb.close <- dist_X_close(bed[,c(1:3)], file.protein.coding );
	idx.distal <- which( tb.close$V7 > 1000)
	idx.proxm <- which( tb.close$V7 <=0 )

cat("<=1K", NROW(idx.proxm), ">1K", NROW(idx.distal), "\n");

	##proximal
	bed.p.ord <- idx.proxm[ order(bed$score[idx.proxm], decreasing=TRUE) ];
	hCountMatrix <- bed.step.bpQuery.bigWig(hMark, center.bed(bed[bed.p.ord ,c(1,5,5)], dist, dist), step=step, abs.value=TRUE)
	hmat.ord <- log(matrix(unlist(hCountMatrix), nrow= NROW(idx.proxm), byrow=TRUE)+1)
	hmat.p <- t(sapply(1:floor(NROW(hmat.ord)/navg), function(x) {colMeans(hmat.ord[((x-1)*navg+1):min(NROW(hmat.ord),(x*navg)),])}))

	##distal
	bed.d.ord <- idx.distal[ order(bed$score[idx.distal], decreasing=TRUE) ];
	hCountMatrix <- bed.step.bpQuery.bigWig(hMark, center.bed(bed[bed.d.ord,c(1,5,5)], dist, dist), step=step, abs.value=TRUE)
	hmat.ord <- log(matrix(unlist(hCountMatrix), nrow= NROW(idx.distal), byrow=TRUE)+1)
	hmat.d <- t(sapply(1:floor(NROW(hmat.ord)/navg), function(x) {colMeans(hmat.ord[((x-1)*navg+1):min(NROW(hmat.ord),(x*navg)),])}))

	png(paste(png.name,".png",sep=""), width=450, height = 800*1.55 )

	nrow.all <- NROW(hmat.p) + NROW(hmat.d);
	lay.heights <- c(NROW(hmat.p)/nrow.all, NROW(hmat.d)/nrow.all, 0.4, 0.15);
	lay.widths  <- c(1, 4, 5)
	layout(matrix(c(1,2,2, 3,4,4, 5,5,5, 6,6,7 ), nrow=4, byrow=T), widths=lay.widths, heights=lay.heights)

	file.tmp <- tempfile(fileext=".bed");
	write.table( bed[bed.p.ord,c(1:3)], file=file.tmp, quote=F, row.names=F, col.names=F, sep="\t");
	tb <- unique(read.table(pipe(paste("zcat ", hPeakFile, "| awk -v OFS='\t' '{print $1,$2,$3 }' - | bedtools intersect -a" , file.tmp, "-b - -loj")))[,c(1:6)]);
	tbovp <- unlist(lapply(1:NROW(tb), function(x){
		if(as.character(tb[x,4]) == ".") return(0)
		else
			return(1);
		#if(tb[x,5] <= tb[x,2] && tb[x,3]<=tb[x,6]) return(1);
		#if(tb[x,3] <= tb[x,5] || tb[x,6]>=tb[x,2]) return(0);
		#if(tb[x,2] < tb[x,5] && tb[x,5] < tb[x,3] && tb[x,3] < tb[x,6] ) return( 1 ) #(tb[x,3]-tb[x,5])/(tb[x,3]-tb[x,2]) );
		#if(tb[x,5] < tb[x,2] && tb[x,2] < tb[x,6] && tb[x,6] < tb[x,3] ) return( 1 ) #(tb[x,6]-tb[x,2])/(tb[x,3]-tb[x,2]) );
		#return(NA);
	}));

	tborg <- bed[bed.p.ord,c(1:3)];
	colnames(tborg) <- c("chr", "start", "end");
	tbovp <- data.frame(tb[,c(1:3)], ratio=tbovp);
	colnames(tbovp) <- c("chr", "start", "end", "ratio");
	tbovp <- sqldf("select chr, start, end, max(ratio) as ratio from tbovp group by chr, start, end")
	tbovp <- sqldf("select tborg.chr, tborg.start, tborg.end, tbovp.ratio from tborg left join tbovp on tborg.chr=tbovp.chr and tborg.start=tbovp.start");

	show(all(tbovp[,c(1:3)]==data.frame(bed[bed.p.ord,c(1:3)], stringsAsFactors=F)))
	ovp.p <- sapply(1:floor(NROW(tbovp)/navg), function(x) {mean(tbovp[((x-1)*navg+1):min(NROW(tbovp),(x*navg)),4])})

	file.tmp <- tempfile(fileext=".bed");
	write.table( bed[bed.d.ord,c(1:3)], file=file.tmp, quote=F, row.names=F, col.names=F, sep="\t");
	tb <- unique(read.table(pipe(paste("zcat ", hPeakFile, "| awk -v OFS='\t' '{print $1,$2,$3 }' - | bedtools intersect -a" , file.tmp, "-b - -loj")))[,c(1:6)]);
	tbovp <- unlist(lapply(1:NROW(tb), function(x){
		if(as.character(tb[x,4]) == ".") return(0)
		else
			return(1);
		#if(tb[x,5] <= tb[x,2] && tb[x,3]<=tb[x,6]) return(1);
		#if(tb[x,3] <= tb[x,5] || tb[x,6]>=tb[x,2]) return(0);
		#if(tb[x,2] < tb[x,5] && tb[x,5] < tb[x,3] && tb[x,3] < tb[x,6] ) return( 1 ) #(tb[x,3]-tb[x,5])/(tb[x,3]-tb[x,2]) );
		#if(tb[x,5] < tb[x,2] && tb[x,2] < tb[x,6] && tb[x,6] < tb[x,3] ) return( 1 ) #(tb[x,6]-tb[x,2])/(tb[x,3]-tb[x,2]) );
		#return(NA);
	}));

	tborg <- bed[bed.d.ord,c(1:3)];
	colnames(tborg) <- c("chr", "start", "end");
	tbovp <- data.frame(tb[,c(1:3)], ratio=tbovp);
	colnames(tbovp) <- c("chr", "start", "end", "ratio");
	tbovp <- sqldf("select chr, start, end, max(ratio) as ratio from tbovp group by chr, start, end")
	tbovp <- sqldf("select tborg.chr, tborg.start, tborg.end, tbovp.ratio from tborg left join tbovp on tborg.chr=tbovp.chr and tborg.start=tbovp.start");

	#show(all(as.character(tbovp[,c(1:3)])==as.character(bed[bed.d.ord,c(1:3)])))
	show(all(tbovp[,c(1:3)]==data.frame(bed[bed.d.ord,c(1:3)], stringsAsFactors=F)))
	ovp.d <- sapply(1:floor(NROW(tbovp)/navg), function(x) {mean(tbovp[((x-1)*navg+1):min(NROW(tbovp),(x*navg)),4])})

	ovp.min <- quantile( c(ovp.p, ovp.d), probs=0.05);
	ovp.max <- quantile( c(ovp.p, ovp.d), probs=0.95);

	##part1:
	par(mar=c(0,0,0,0), plt=c(0.2, 0.8,0.01, 0.99 ));
	plot(NA,NA, type="n", xlab="", ylab="", xlim=c(0,1), ylim=c(0, NROW(ovp.p)), xaxs="i", yaxs="i", xaxt = "n", yaxt = "n", bty="n" );
	for(i in 1:NROW(ovp.p))
	{
		#cl.base <- ovp.p[NROW(ovp.p)+1-i];
		cl.base <- (ovp.p[NROW(ovp.p)+1-i] - ovp.min)/(ovp.max- ovp.min);
		if(cl.base<0) cl.base <- 0;
		if(cl.base>1) cl.base <- 1;
		clCor <- rgb(1-cl.base, 1-cl.base, 1-cl.base)
		rect(0, i-1, 1, i, col = clCor, border = clCor)
	}

	##part2:
	##png(paste(png.name,".proxmial.png",sep=""), width=400, height = round(800/NROW(hmat.all)*NROW(hmat.p)) )     # width and height are in pixels
	gt <- pheatmap( hmat.p, cluster_rows = FALSE, cluster_cols = FALSE, col= hmcols, breaks = bk, legend=FALSE, show_rownames=FALSE, show_colnames=FALSE, silent=T )

	par(mar=c(0,0,0,0), plt=c(0.2, 0.8,0.2, 0.8 ));
    plot(NA,NA, type="n", xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i", xaxt = "n", yaxt = "n");
    ##grid.newpage()
	pushViewport(viewport(layout = grid.layout(4, 3, widths=lay.widths, heights=lay.heights) ))
	gt$gtable$vp <- viewport(layout.pos.row = 1, layout.pos.col = c(2,3))
    grid.draw(gt$gtable)
	popViewport()

	##part3:
	par(mar=c(0,0,0,0), plt=c(0.2, 0.8,0.01, 0.99 ));
	plot(NA,NA, type="n", xlab="", ylab="", xlim=c(0,1), ylim=c(0, NROW(ovp.d)), xaxs="i", yaxs="i", xaxt = "n", yaxt="n", bty="n" );
	for(i in 1:NROW(ovp.d))
	{
		cl.base <- (ovp.d[NROW(ovp.d)+1-i] - ovp.min)/(ovp.max- ovp.min);
		if(cl.base<0) cl.base <- 0;
		if(cl.base>1) cl.base <- 1;
		clCor <- rgb( 1-cl.base, 1-cl.base, 1-cl.base)
		rect(0, i-1, 1, i, col = clCor, border = clCor);
	}

	##part4:
	##png(paste(png.name,".distal.png",sep=""), width=400, height = round(800/NROW(hmat.all)*NROW(hmat.d)))     # width and height are in pixels
	gt <- pheatmap( hmat.d, cluster_rows = FALSE, cluster_cols = FALSE, col= hmcols, breaks = bk, legend=FALSE, show_rownames=FALSE, show_colnames=FALSE, silent=T )
	par(mar=c(0,0,0,0), plt=c(0.2, 0.8,0.2, 0.8 ));
    plot(NA,NA, type="n", xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i", xaxt = "n", yaxt = "n");
    ##grid.newpage()
	pushViewport(viewport(layout = grid.layout(4, 3, widths=lay.widths, heights=lay.heights)))
	gt$gtable$vp <- viewport(layout.pos.row = 2, layout.pos.col =  c(2,3))
	grid.draw(gt$gtable)
	popViewport()

	##part 5:
	## draw mat-plot for overlapping peak and non-overlapped peak
	file.tmp <- tempfile(fileext=".bed");
	write.table( bed[,c(1:3)], file=file.tmp, quote=F, row.names=F, col.names=F, sep="\t");
	tb <- unique(read.table(pipe(paste("zcat ", hPeakFile, " | bedtools intersect -a" , file.tmp, "-b - -loj")))[,c(1:4)]);
	#no-overlap
	idx.nonoverlap <- which( unlist(apply(tb, 1, function(x){return(as.character(x[4])==".")})) );

	hCountMatrix <- bed.step.bpQuery.bigWig(hMark, center.bed(bed[,c(1,5,5)], dist, dist), step=step, abs.value=TRUE)
	hmat <- matrix(unlist(hCountMatrix), nrow= NROW(bed), byrow=TRUE);
	mat.p <- colMeans(hmat[idx.nonoverlap,]);
	mat.m <- colMeans(hmat[-idx.nonoverlap,]);

	par(mar=c(10,0,0,0), plt=c(0.125, 0.99, 0.2, 0.90), mgp=c(5,2,0) );
	##png( paste(png.name, ".metaplot.png",sep=""), width=400, height=300  )
	y.max <- max(c(mat.p, mat.m));
	plot(NA, NA, type="n", col="gray", xlim=c(1,NROW(mat.p)), ylim=c(0, y.max)/1000,  xlab="Distance (Kbp) ", ylab="", cex.axis=3, cex.lab=3,  xaxt = "n", yaxt="n")
	## overlap peak
	lines(1:NROW(mat.m), mat.m/1000, col="darkgreen", lwd=2);
	## non-overlap peak
	lines(1:NROW(mat.p), mat.p/1000, col="blue", lwd=2);
	axis(1, c(0, 100, 200, 300, 400), c(-dist/1000, -dist/2000, 0, dist/2000, dist/1000), cex.axis=3, cex=3 )
	if( y.max>10 )
		axis(2, c(0, y.max/2000, y.max/1000), c(0, round(y.max/2), round(y.max)), cex.axis=3, cex=3 )
	else if( y.max>=1 )
		axis(2, c(0, y.max/2000, y.max/1000), c(0, round(y.max/2, digits=1), round(y.max,digits=1)), cex.axis=3, cex=3 )
	else
		axis(2, c(0, y.max/2000, y.max/1000), c(0, signif(y.max/2,digits=2), signif(y.max,digits=2)), cex.axis=3, cex=3 )

	##part 6:
	## draw colorScale for peak track
	par(mar=c(10,0,0,0), plt=c(0.15, 0.85, 0.4, 0.8 ));
	plot(NA,NA, type="n", xlim=c(0, 1), ylim=c(0,1), xlab="", ylab="", bty="n", xaxt="n", yaxt= "n" );
	for(i in seq(0, 0.95, 0.05) )
	{
		clCor <- rgb( 1-i, 1-i, 1-i)
		rect(i,0, i+0.05, 1, col=clCor, border=clCor);
	}
	axis(1, c(0, 0.95), c(signif(ovp.min, digits=1), signif(ovp.max, digits=1) ), cex.axis=3, cex=3, tick=FALSE );

	##part 7:
	## draw colorScale for Heatmap
    if(is.null(bk))
       bk <- seq(min(hmat), max(hmat), 0.01)
    if(is.null(hmcols))
       hmcols <- rev(colorRampPalette(brewer.pal(9,"RdBu"))(length(bk)-1))

	par(mar=c(10,0,0,0), plt=c(0.15, 1.00, 0.4, 0.8 ));
	plot(NA,NA, type="n", xlim=c(1, NROW(bk)*0.85/0.7), ylim=c(0,1), xlab="", ylab="", bty="n", xaxt="n", yaxt= "n" );
	for(i in 1:NROW(bk)) rect(i,0, i+1, 1, col=hmcols[i], border=hmcols[i]);
	axis(1, c(1, NROW(bk)/2, NROW(bk)), round(exp(c(bk[1], bk[round(NROW(bk)/2)], bk[NROW(bk)])),0), cex.axis=3, cex=3, tick=FALSE );

	dev.off();

	return(list(breaks= bk, hmcols=hmcols))
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

if(0)
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

	file.dnase.bw = "/fs/cbsudanko/storage/data/hg19/k562/dnase/wgEncodeUwDgfK562Raw.bigWig" #"wgEncodeOpenChromDnaseK562SigV2.bigWig";
	#file.dnase.peak = "/fs/cbsudanko/storage/data/hg19/k562/dnase/GSM646567_hg19_wgEncodeUwDgfK562Pk.narrowPeak.txt.gz" #"wgEncodeOpenChromDnaseK562PkV2.narrowPeak.gz";
	file.dnase.peak = "./GSM646567_hg19_wgEncodeUwDgfK562Pk.macs2.narrowPeak.gz";

	file.grocap <- "../k562/hg19.k562.grocap.pair.bed";
	file.protein.coding <- "./gencode.v19.protein.coding.bed"

	#file.plus.bw <- "../k562/6045_7157_27176_HNHKJBGXX_K562_0min_celastrol10uM_rep2_GB_CAGATC_R1_plus.primary.bw"
	#file.minus.bw <- "../k562/6045_7157_27176_HNHKJBGXX_K562_0min_celastrol10uM_rep2_GB_CAGATC_R1_minus.primary.bw"
	file.plus.bw <- "../k562/K562_unt.sort.bed.gz_plus.bw"
	file.minus.bw <- "../k562/K562_unt.sort.bed.gz_minus.bw"

	file.prefix <- "dregG1";
	dregX <- read.table("../new-rf-201803/G1/G1.dREG.peak.full.bed.gz");
	dregX <- dregX[dregX$V5<=0.05,];
	dregX <- dregX[,-5];
	colnames(dregX) <- c("chr", "start", "end", "score", "peak")
	dregX <- dregX[grep("_|chrY", dregX[,1], invert=TRUE),]
	dregX <- cbind(dregX, width=dregX$end-dregX$start);

	sup <- writeHeatmap( dregX, file.plus.bw, file.minus.bw, file.H3K27ac.bw, file.H3K27ac.peak, paste(file.prefix, "-H3K27ac",sep="") )
	##writeColorScale(sup$breaks, cols=sup$hmcols, paste(file.prefix, "-H3K27ac",sep="") )

	sup <- writeHeatmap( dregX, file.plus.bw, file.minus.bw, file.dnase.bw,   file.dnase.peak,   paste(file.prefix, "-DNase",sep="") )
	##writeColorScale(sup$breaks, cols=sup$hmcols, paste(file.prefix, "-DNase",sep="") )

	sup <- writeHeatmap( dregX, file.plus.bw, file.minus.bw, file.H3K4me1.bw, file.H3K4me1.peak, paste(file.prefix, "-H3K4me1",sep="") )
	##writeColorScale(sup$breaks, cols=sup$hmcols, paste(file.prefix, "-H3K4me1",sep="") )

	sup <- writeHeatmap( dregX, file.plus.bw, file.minus.bw, file.H3K4me3.bw, file.H3K4me3.peak, paste(file.prefix, "-H3K4me3",sep="") )
	##writeColorScale(sup$breaks, cols=sup$hmcols, paste(file.prefix, "-H3K4me3",sep="") )

	sup <- writeHeatmap( dregX, file.plus.bw, file.minus.bw, file.H3K27me3.bw,file.H3K27me3.peak, paste(file.prefix, "-H3K27me3",sep="") )
	##writeColorScale(sup$breaks, cols=sup$hmcols, paste(file.prefix, "-H3K27me3",sep="") )

	##sup <- writeHeatmap( dregX, file.plus.bw, file.minus.bw, file.H3K9me3.bw, file.H3K9me3.peak, paste(file.prefix, "-H3K9me3",sep="") )
	##writeColorScale(sup$breaks, cols=sup$hmcols, paste(file.prefix, "-H3K9me3",sep="") )

	##sup <- writeHeatmap( dregX, file.MNase.bw, paste(file.prefix, "-MNase",sep=""), hm_order= ord, path=MNase.path)
	##writeColorScale(sup$breaks, cols=NULL, paste(file.prefix, "-MNase",sep="") )
}

writeMatplot<- function(bed, file.plus.bw, file.minus.bw, hMarkFile, hPeakFile, title, subs= NULL, breaks= NULL,
							cols= NULL, dist= 5000, step=25)
{
	## Load mark.
	hMark <- load.bigWig(hMarkFile)

	## draw mat-plot for overlapping peak and non-overlapped peak
	file.tmp <- tempfile(fileext=".bed");
	write.table( bed[,c(1:3)], file=file.tmp, quote=F, row.names=F, col.names=F, sep="\t");
	tb <- unique(read.table(pipe(paste("zcat ", hPeakFile, " | bedtools intersect -a" , file.tmp, "-b - -loj")))[,c(1:4)]);
	idx.overlap <- which( unlist(apply(tb, 1, function(x){return(as.character(x[4])!=".")})) );

	hCountMatrix <- bed.step.bpQuery.bigWig(hMark, center.bed(bed[,c(1,5,5)], dist, dist), step=step, abs.value=TRUE)
	hmat <- matrix(unlist(hCountMatrix), nrow= NROW(bed), byrow=TRUE);
	mat.p <- colMeans(hmat[idx.overlap,]);
	mat.m <- colMeans(hmat[-idx.overlap,]);

	par( mar=c(5,4,2,2), plt=c(0.15, 0.99, 0.2, 0.8), mgp=c(1.5,0.4,0) );

	#y.max <- max(c( mat.m));
	y.max <- 150;
	plot(NA, NA, type="n", col="gray", xlim=c(1,NROW(mat.p)), ylim=c(0, y.max)/1000,  xlab="Distance (Kbp) ", ylab="", cex.axis=1, cex.lab=1,  xaxt = "n", yaxt="n", main=title)
	lines(1:NROW(mat.m), mat.m/1000, col="darkgreen", lwd=1);
	#lines(1:NROW(mat.p), mat.p/1000, col="blue", lwd=1);
	axis(1, c(0, 100, 200, 300, 400), c(-dist/1000, -dist/2000, 0, dist/2000, dist/1000), cex.axis=1, cex=1 )
	if( y.max>10 )
		axis(2, c(0, y.max/2000, y.max/1000), c(0, round(y.max/2), round(y.max)), cex.axis=1, cex=1 )
	else if( y.max>=1 )
		axis(2, c(0, y.max/2000, y.max/1000), c(0, round(y.max/2, digits=1), round(y.max,digits=1)), cex.axis=1, cex=1 )
	else
		axis(2, c(0, y.max/2000, y.max/1000), c(0, signif(y.max/2,digits=2), signif(y.max,digits=2)), cex.axis=1, cex=1 )


	return()
}

if(0)
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

	file.dnase.bw = "/fs/cbsudanko/storage/data/hg19/k562/dnase/wgEncodeUwDgfK562Raw.bigWig" #"wgEncodeOpenChromDnaseK562SigV2.bigWig";
	#file.dnase.peak = "/fs/cbsudanko/storage/data/hg19/k562/dnase/GSM646567_hg19_wgEncodeUwDgfK562Pk.narrowPeak.txt.gz" #"wgEncodeOpenChromDnaseK562PkV2.narrowPeak.gz";
	file.dnase.peak = "./GSM646567_hg19_wgEncodeUwDgfK562Pk.macs2.narrowPeak.gz";

	file.grocap <- "../k562/hg19.k562.grocap.pair.bed";
	file.protein.coding <- "./gencode.v19.protein.coding.bed"

	#file.plus.bw <- "../k562/6045_7157_27176_HNHKJBGXX_K562_0min_celastrol10uM_rep2_GB_CAGATC_R1_plus.primary.bw"
	#file.minus.bw <- "../k562/6045_7157_27176_HNHKJBGXX_K562_0min_celastrol10uM_rep2_GB_CAGATC_R1_minus.primary.bw"
	file.plus.bw <- "../k562/K562_unt.sort.bed.gz_plus.bw"
	file.minus.bw <- "../k562/K562_unt.sort.bed.gz_minus.bw"

	dregX <- read.table("../new-rf-201803/G1/G1.dREG.peak.full.bed.gz");

	pdf("fig-S10.pdf",width=4.5, height=2);

	dregX <- dregX[dregX$V5<=0.05,];
	dregX <- dregX[, -5 ];
	colnames(dregX) <- c("chr", "start", "end", "score", "peak")
	dregX <- dregX[grep("_|chrY", dregX[,1], invert=TRUE),]
	dregX <- cbind(dregX, width=dregX$end-dregX$start);

	layout(matrix(c(1, 2, 3 ), nrow=1, byrow=T), widths=rep(1/3,3), heights=1)

	#sup <- writeMatplot( dregX, file.plus.bw, file.minus.bw, file.dnase.bw,   file.dnase.peak, "dNase")
	sup <- writeMatplot( dregX, file.plus.bw, file.minus.bw, file.H3K27ac.bw, file.H3K27ac.peak, "H3k27ac")
	sup <- writeMatplot( dregX, file.plus.bw, file.minus.bw, file.H3K4me1.bw, file.H3K4me1.peak, "H3k4me1")
	sup <- writeMatplot( dregX, file.plus.bw, file.minus.bw, file.H3K4me3.bw, file.H3K4me3.peak, "H3k4me3")

	dev.off();
}

if(0)
{
	pdf("G1-peak-dist.pdf");
	slices <- c(15652, 74115-43127-15652, 43127 )
	lbls <- c("Proximal", "<=1000", "Distal")
	pie(slices, labels = lbls, main="Distance to TSS");
	dev.off()
}
