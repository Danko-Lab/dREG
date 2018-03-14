library(bigWig);

#/fs/cbsudanko/storage/data/hg19/k562/histones/wgEncodeBroadHistoneK562H3k27acStdAln.bed.gz
TPR.H3k27ac <- c(0.783, 0.526, 0.733, 0.635, 0.687, 0.582, 0.581, 0.718)
#../k562/hg19.k562.grocap.pair.bed
TPR.grocap <-  c(0.945, 0.782, 0.923, 0.840, 0.894, 0.812, 0.813, 0.961)

#/fs/cbsudanko/storage/data/hg19/k562/histones/wgEncodeBroadHistoneK562H3k27acStdAln.bed.gz
#TPR.H3k27ac <- c(0.8027,0.5486,0.7470,0.6707,0.7035,0.5959,0.6006,0.8346)
##../k562/hg19.k562.grocap.pair.bed
#TPR.grocap <- c(0.963, 0.820, 0.940, 0.837, 0.916, 0.840, 0.843, 0.970)

dataname <- c(paste("K562-G", 1:7, sep=""), "GM12878")
cols <- c("#a70ca0", "#86db63", "#5070fb", "#e27800", "#f7aaff", "#017657", "#f92367", "#8e2f59" );


file.bw.G1 <- c( "../k562/K562_unt.sort.bed.gz_plus.bw", "../k562/K562_unt.sort.bed.gz_minus.bw");
file.bw.G2 <- c( "../k562/groseq_plus.bigWig", "../k562/groseq_minus.bigWig");
file.bw.G3 <- c( "../k562/K562_Nuc_NoRNase_plus.bw", "../k562/K562_Nuc_NoRNase_minus.bw");
file.bw.G4 <- c( "../k562/K562_Nuc_RNase_plus.bw", "../k562/K562_Nuc_RNase_minus.bw");
file.bw.G5 <- c( "../k562/K562_FC_NHS_BRs_normalized_pl.bigWig", "../k562/K562_FC_NHS_BRs_normalized_mn.bigWig");
file.bw.G6 <- c( "../k562/6045_7157_27170_HNHKJBGXX_K562_0min_celastrol10uM_rep1_GB_ATCACG_R1_plus.primary.bw", "../k562/6045_7157_27170_HNHKJBGXX_K562_0min_celastrol10uM_rep1_GB_ATCACG_R1_minus.primary.bw");
file.bw.G7 <- c( "../k562/6045_7157_27176_HNHKJBGXX_K562_0min_celastrol10uM_rep2_GB_CAGATC_R1_plus.primary.bw", "../k562/6045_7157_27176_HNHKJBGXX_K562_0min_celastrol10uM_rep2_GB_CAGATC_R1_minus.primary.bw");
file.bw.GM <- c( "../GM12878/groseq_plus.bigWig", "../GM12878/groseq_minus.bigWig");

file.K562.H3K27ac.peak  <- "/fs/cbsudanko/storage/data/hg19/k562/histones/wgEncodeBroadHistoneK562H3k27acStdAln.bed.gz"
#file.K562.DHS.peak  <- "../k562/K562.bed";
file.K562.DHS.peak  <- "./k562.merged.broad.peak.bed";
file.grocap="../k562/hg19.k562.new_hmm2b.post2.bed"

make_big_bed<-function(file.nonneg )
{
	system( paste("cat ../k562/wgEncodeOpenChromDnaseK562PkV2.narrowPeak ../k562/GSM646567_hg19_wgEncodeUwDgfK562Pk.narrowPeak.txt ../k562/GSM646567_hg19_wgEncodeUwDgfK562Pk.macs2.narrowPeak | awk -v OFS='\\t' '{print $1,$2,$3}' - > ", file.nonneg ), intern=TRUE );

	nonneg_bed <- unique(read.table( pipe(paste("cat ", file.nonneg, file.grocap, " | sort-bed ", file.nonneg," | bedtools merge -i - -d 100 ", sep=" " ) ) )[,c(1:3)]);
	nonneg_bed[,2] <- nonneg_bed[,2] - 100
	idx.mis <- which(nonneg_bed[,2]<0);
	if(length(idx.mis)>0) nonneg_bed[idx.mis,2] <- 0;
	nonneg_bed[,3] <- nonneg_bed[,3] + 100
	write.table( nonneg_bed, file=file.nonneg, quote=F, row.name=F, col.names=F, sep="\t" );
}

make_big_bed( file.K562.DHS.peak ); 



file.GM.H3K27ac.peak  <- "/fs/cbsudanko/storage/data/hg19/gm12878/histones/wgEncodeBroadHistoneGm12878H3k27acStdPk.broadPeak.gz"
file.GM12878.DHS.peak <- "/fs/cbsudanko/storage/data/hg19/gm12878/dnase/uw.merge.narrowPeak.bed"

file.dREG.G1 <- "/workdir/zw355/proj/prj10-dreg/new-rf-201803/G1/G1.dREG.peak.score.bed.gz";
file.dREG.G2 <- "/workdir/zw355/proj/prj10-dreg/new-rf-201803/G2/G2.dREG.peak.score.bed.gz";
file.dREG.G3 <- "/workdir/zw355/proj/prj10-dreg/new-rf-201803/G3/G3.dREG.peak.score.bed.gz";
file.dREG.G4 <- "/workdir/zw355/proj/prj10-dreg/new-rf-201803/G4/G4.dREG.peak.score.bed.gz";
file.dREG.G5 <- "/workdir/zw355/proj/prj10-dreg/new-rf-201803/G5/G5.dREG.peak.score.bed.gz";
file.dREG.G6 <- "/workdir/zw355/proj/prj10-dreg/new-rf-201803/G6/G6.dREG.peak.score.bed.gz";
file.dREG.G7 <- "/workdir/zw355/proj/prj10-dreg/new-rf-201803/G7/G7.dREG.peak.score.bed.gz";
file.dREG.GM <- "/workdir/zw355/proj/prj10-dreg/new-rf-201803/GM/GM.dREG.peak.score.bed.gz";
file.dREG <-c( file.dREG.G1, file.dREG.G2, file.dREG.G3, file.dREG.G4, file.dREG.G5, file.dREG.G6, file.dREG.G7, file.dREG.GM);

file.bws <- list(file.bw.G1, file.bw.G2, file.bw.G3, file.bw.G4, file.bw.G5, file.bw.G6, file.bw.G7, file.bw.GM);
coverage <- basesCovered <- datasize <- c();
for(i in 1:length(file.bws))
{
	bw1 <- load.bigWig(file.bws[[i]][1]);
	bw2 <- load.bigWig(file.bws[[i]][2]);
	coverage <-c(coverage, abs(bw1$mean*bw1$basesCovered) + abs(bw2$mean*bw2$basesCovered));
	basesCovered<-c(basesCovered, bw1$basesCovered+bw2$basesCovered );
	datasize<-c(datasize, bw1$primaryDataSize+bw2$primaryDataSize );
	unload.bigWig(bw1);
	unload.bigWig(bw2);
}

n.overlap <- n.novel<-c();
for(i in 1:length(file.dREG ) )
{
	if(i==8)
	{
		tb <- unique(read.table(pipe(paste("zcat ", file.GM.H3K27ac.peak , " | bedtools intersect -a" , file.dREG[i], "-b - -loj")))[,c(1:5)])
		idx.novel <- which( unlist(apply(tb, 1, function(x){return(as.character(x[5])==".")})) );
		tb.nov <- tb [ idx.novel, 1:4];
		
		
		file1= tempfile();
		write.table(tb.nov, file=file1, quote=F, row.names=F, col.names=F, sep="\t");

		tb <- unique(read.table(pipe(paste("cat ", file.GM12878.DHS.peak, " | bedtools intersect -a" , file1, "-b - -loj")))[,c(1:5)]);
		idx.novel <- which( unlist(apply(tb, 1, function(x){return(as.character(x[5])==".")})) );
		tb.nov <- tb[ idx.novel, 1:4];
		
		unlink(file1);

	}
	else
	{
		tb <- unique(read.table(pipe(paste("zcat ", file.K562.H3K27ac.peak , " | bedtools intersect -a" , file.dREG[i], "-b - -loj")))[,c(1:5)]);
		idx.novel <- which( unlist(apply(tb, 1, function(x){return(as.character(x[5])==".")})) );
		tb.nov <- tb [ idx.novel, 1:4]
		
		file1= tempfile();
		write.table(tb.nov, file=file1, quote=F, row.names=F, col.names=F, sep="\t");

		tb <- unique(read.table(pipe(paste("cat ", file.K562.DHS.peak, " | bedtools intersect -a" , file1, "-b - -loj")))[,c(1:5)]);
		idx.novel <- which( unlist(apply(tb, 1, function(x){return(as.character(x[5])==".")})) );
		tb.nov <- tb[ idx.novel, 1:4];
		
		unlink(file1);
	}
	

	n.novel <- c( n.novel, NROW(tb.nov) );
	n.overlap <- c(n.overlap, NROW(tb)-NROW(tb.nov) );
}


TPR.H3k27ac <- TPR.H3k27ac[-4];
TPR.grocap <- TPR.grocap[-4];
basesCovered <- basesCovered[-4];
coverage <- coverage[-4];
datasize <- datasize[-4];
dataname <- dataname[-4]
cols <- cols[-4]

n.novel <- n.novel[-4]
n.overlap <- n.overlap[-4]

cols <- c("#86db63", "#86db63", "#86db63", "#86db63", "#86db63", "#f92367", "#5070fb" );

show(coverage);
#log.coverage <- log(coverage)
coverage <- coverage/1000/1000

show(n.novel);

pdf("fig-2B.pdf")

plot( coverage, n.novel, main="", xlab="Number of mapped reads (M)", ylab="#Novel elements", ylim=c(0,max(n.novel)), pch=17, col=cols, cex=2)
r <- summary(lm( n.novel~ log(coverage)));
r.b <- r$coefficients[1,1]
r.k <- r$coefficients[2,1]
#r.x <- c( seq( min(log.coverage), max(log.coverage), 0.1 ),  max(log.coverage))
r.x <- seq( 1, 500, 1 );
lines( r.x, r.b+r.k*log(r.x), lty=22, lwd=1, col="black");

#axis(1, at=seq(16.5, 20, 0.7), labels=round( exp(seq(16.5, 20, 0.7))/1000/1000, 1), cex=2);

#points(basesCovered/1000000, n.overlap, pch=19, col=cols, cex=2)
#r <- summary(lm(n.overlap~basesCovered));
#r.b <- r$coefficients[1,1]
#r.k <- r$coefficients[2,1]
#r.x <- c( seq( min(basesCovered), max(basesCovered), 1000000 ),  max(basesCovered))
#lines( r.x/1000000, r.b+r.k*r.x, lty=22, lwd=1, col="black");

legend.name <- c( "Training", "K562 holdout", "GM12878 holdout");
legend.cols  <- c( "#86db63", "#f92367", "#5070fb")
legend.pch  <- c( 17, 17, 17)
legend("bottomright", legend.name, text.col=legend.cols, col=legend.cols, pch=legend.pch, cex=1, ncol=1, horiz=F);

dev.off();

