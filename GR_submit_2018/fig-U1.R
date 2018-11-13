library(bigWig)

file.grocap <- "../k562/hg19.k562.grocap.pair.bed";

dist_tss_center <- function( pred.bed )
{
	options("scipen"=100, "digits"=4);

	tmp.bed <- tempfile(fileext=".bed" )
	write.table( pred.bed[,c(1:3)], file=tmp.bed, row.names=F, col.names=F, quote=F, sep="\t");
	system( paste("sort-bed ",  tmp.bed, " > ", tmp.bed, ".sorted", sep="") );

	tb.grocap <- read.table( file.grocap, header=F);
	#tb.grocap[,2] <- round((tb.grocap[,2] + tb.grocap[,3])/2)
	#tb.grocap[,3] <- tb.grocap[,2] + 1;
	tmp2.bed <- tempfile(fileext=".bed")
	write.table(tb.grocap, file=tmp2.bed, row.names=F, col.names=F, quote=F, sep="\t");

	tb.close <-  read.table( file = pipe(paste("sort-bed ",  tmp2.bed, " | bedtools closest -d -a ", tmp.bed, ".sorted -b - -t first", sep="") ) );
	return(tb.close);
}


draw_peak_dist <- function( title, file.dreg.peak, file.olddreg.peak, file.Tfit.peak, break.width=50, plot.width=1000)
{
	dregP <- read.table(file.dreg.peak, header=F);
	dregX <- dregP[,c(1, 8, 8)];
	dregX[,3] <- dregX[,3] + 1;
	tb.close1 <- dist_tss_center(dregX);
	tb.dist <- apply(tb.close1[,c(2,5,6)], 1, function(x){ d<-0; if(x[1]<x[2]) d <- x[1]- x[2]; if(x[1]>x[3]) d <- x[1]- x[3]; d });
	x1 <- hist(tb.dist, breaks=max(tb.dist)/break.width, plot=FALSE);

	dregX <- read.table(file.olddreg.peak, header=F );
	dregX[,2] <- round((dregX[,2] + dregX[,3])/2)
	dregX[,3] <- dregX[,2] + 1;
	tb.close4 <- dist_tss_center(dregX);
	tb.dist <- apply(tb.close4[,c(2,5,6)], 1, function(x){ d<-0; if(x[1]<x[2]) d <- x[1]- x[2]; if(x[1]>x[3]) d <- x[1]- x[3]; d });
	x2 <- hist(tb.dist, breaks=max(tb.dist)/break.width, plot=FALSE);

	dregX <- read.table(file.Tfit.peak, header=F );
	dregX[,2] <- round((dregX[,2] + dregX[,3])/2)
	dregX[,3] <- dregX[,2] + 1;
	tb.close7 <- dist_tss_center(dregX);
	tb.dist <- apply(tb.close7[,c(2,5,6)], 1, function(x){ d<-0; if(x[1]<x[2]) d <- x[1]- x[2]; if(x[1]>x[3]) d <- x[1]- x[3]; d });
	x3 <- hist(tb.dist, breaks=max(tb.dist)/break.width, plot=FALSE);

	plot(x1$breaks[-1], x1$counts, type="l",  col="red", xlim=c(-1,1)*plot.width, main=title, xlab="Peak distance", ylab="Peak Counts", lwd=2);
	lines(x2$breaks[-1], x2$counts, col="black", lwd=2);
	lines(x3$breaks[-1], x3$counts, col="darkblue", lwd=2);
	legend("topright", c("new dREG", "old dREG", "Tfit divergent"),
		text.col=c("red", "black", "darkblue"), lwd=2, col=c("red", "black", "darkblue"))
}

if(0)
{
pdf("fig-U1-dist.pdf");
draw_peak_dist( "K562 holdout", "../new-model-201709/G7/out.dREG.peak.full.bed.gz", "../old-model/test-G7/out.dREG.peak.gz", "~/temp/Gx/6045_rep2.Tfit.divergent.bed.gz", break.width=20, plot.width=1000);
draw_peak_dist( "GM12878 holdout", "../new-model-201709/GM12878/out.dREG.peak.full.bed.gz", "../old-model/test-GM12878/out.dREG.peak.gz", "../Tfit/gm12878_groseq-2_divergent_classifications.bed", break.width=20, plot.width=1000);
dev.off();
}

draw_peak_reads <- function( title, file.plus.bw, file.minus.bw, file.dreg.peak, file.olddreg.peak, file.Tfit.peak, break.width=50, plot.width=1000)
{
	bw.plus<-load.bigWig(file.plus.bw)
	bw.minus<-load.bigWig(file.minus.bw)
	scaling_func<-function(x){x/ break.width}

	dregP <- read.table(file.dreg.peak, header=F);
	dregX <- dregP[,c(1, 8, 8)];
	dregX[,2] <- dregX[,2] - plot.width;
	dregX[,3] <- dregX[,3] + plot.width;

	p1 <- metaprofile.bigWig( dregX, bw.plus, bw.plus, step=break.width, matrix.op= scaling_func)
	m1 <- metaprofile.bigWig( dregX, bw.minus, bw.minus, step=break.width, matrix.op= scaling_func)

	dregX <- read.table(file.olddreg.peak, header=F );
	dregX[,2] <- round((dregX[,2] + dregX[,3])/2)
	dregX[,3] <- dregX[,2];
	dregX[,2] <- dregX[,2] - plot.width;
	dregX[,3] <- dregX[,3] + plot.width;

	p2 <- metaprofile.bigWig( dregX, bw.plus, bw.plus, step=break.width, matrix.op= scaling_func)
	m2 <- metaprofile.bigWig( dregX, bw.minus, bw.minus, step=break.width, matrix.op= scaling_func)

	dregX <- read.table(file.Tfit.peak, header=F );
	dregX[,2] <- round((dregX[,2] + dregX[,3])/2)
	dregX[,3] <- dregX[,2];
	dregX[,2] <- dregX[,2] - plot.width;
	dregX[,3] <- dregX[,3] + plot.width;

	p3 <- metaprofile.bigWig( dregX, bw.plus, bw.plus, step=break.width, matrix.op= scaling_func)
	m3 <- metaprofile.bigWig( dregX, bw.minus, bw.minus, step=break.width, matrix.op= scaling_func)

	plot(1,1, type="n",  col="red", xlim=c(1,NROW(p1$middle)), ylim=c(min(c(m1$middle, m2$middile, m3$middle)), max(c(p1$middle, p2$middile, p3$middle))), main=title, xlab="Peak distance", ylab="reads per 50bp");

	for(i in 1:NROW(p1$middle))
	{
		lines(c(i,i), c(p1$top[i], p1$bottom[i]), col= adjustcolor( "red", alpha.f = 0.2), lwd=0.5 )
		lines(c(i,i), c(m1$top[i], m1$bottom[i]), col= adjustcolor( "red", alpha.f = 0.2), lwd=0.5 )
	}
	lines(1:NROW(p1$middle), p1$middle, col="red", lwd=2);
	lines(1:NROW(m1$middle), m1$middle, col="red", lwd=2);
	

	for(i in 1:NROW(p2$middle))
	{
		lines(c(i,i), c(p2$top[i], p2$bottom[i]), col= adjustcolor( "black", alpha.f = 0.2), lwd=0.5 )
		lines(c(i,i), c(m2$top[i], m2$bottom[i]), col= adjustcolor( "black", alpha.f = 0.2), lwd=0.5 )
	}
	lines(1:NROW(p2$middle), p2$middle, col="black", lwd=2);
	lines(1:NROW(m2$middle), m2$middle, col="black", lwd=2);

	for(i in 1:NROW(p3$middle))
	{
		lines(c(i,i), c(p3$top[i], p3$bottom[i]), col= adjustcolor( "darkblue", alpha.f = 0.2), lwd=0.5 )
		lines(c(i,i), c(m3$top[i], m3$bottom[i]), col= adjustcolor( "darkblue", alpha.f = 0.2), lwd=0.5 )
	}
	lines(1:NROW(p3$middle), p3$middle, col="darkblue", lwd=2);
	lines(1:NROW(m3$middle), m3$middle, col="darkblue", lwd=2);

	legend("topright", c("new dREG", "old dREG", "Tfit divergent"),
		text.col=c("red", "black", "darkblue"), lwd=2, col=c("red", "black", "darkblue"))
}


file.plus.bw <- "../k562/6045_7157_27176_HNHKJBGXX_K562_0min_celastrol10uM_rep2_GB_CAGATC_R1_plus.primary.bw"
file.minus.bw <- "../k562/6045_7157_27176_HNHKJBGXX_K562_0min_celastrol10uM_rep2_GB_CAGATC_R1_minus.primary.bw"

if(1)
{
pdf("fig-U1-dist.pdf");
draw_peak_reads( "K562 holdout", file.plus.bw, file.minus.bw, "../new-model-201709/G7/out.dREG.peak.full.bed.gz", "../old-model/test-G7/out.dREG.peak.gz", "~/temp/Gx/6045_rep2.Tfit.divergent.bed.gz", break.width=50, plot.width=2000);
#draw_peak_reads( "GM12878 holdout", "../new-model-201709/GM12878/out.dREG.peak.full.bed.gz", "../old-model/test-GM12878/out.dREG.peak.gz", "../Tfit/gm12878_groseq-2_divergent_classifications.bed", break.width=20, plot.width=1000);
dev.off();
}
