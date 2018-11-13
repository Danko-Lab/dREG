library(bigWig);

file.black.bed <- "/fs/cbsudanko/storage/data/hg19/all/encode_blacklist/hg19.blacklist.bed";

dataname <- c(paste("K562-G", 1:8, sep=""), "GM12878", "HCT116", "CD4", "MCF7", "HELA" )

file.bw.G1 <- c( "../k562/K562_unt.sort.bed.gz_plus.bw", "../k562/K562_unt.sort.bed.gz_minus.bw");
file.bw.G2 <- c( "../k562/groseq_plus.bigWig", "../k562/groseq_minus.bigWig");
file.bw.G3 <- c( "../k562/K562_Nuc_NoRNase_plus.bw", "../k562/K562_Nuc_NoRNase_minus.bw");
file.bw.G4 <- c( "../k562/K562_Nuc_RNase_plus.bw", "../k562/K562_Nuc_RNase_minus.bw");
file.bw.G5 <- c( "../k562/K562_FC_NHS_BRs_normalized_pl.bigWig", "../k562/K562_FC_NHS_BRs_normalized_mn.bigWig");
file.bw.G6 <- c( "../k562/6045_7157_27170_HNHKJBGXX_K562_0min_celastrol10uM_rep1_GB_ATCACG_R1_plus.primary.bw", "../k562/6045_7157_27170_HNHKJBGXX_K562_0min_celastrol10uM_rep1_GB_ATCACG_R1_minus.primary.bw");
file.bw.G7 <- c( "../k562/6045_7157_27176_HNHKJBGXX_K562_0min_celastrol10uM_rep2_GB_CAGATC_R1_plus.primary.bw", "../k562/6045_7157_27176_HNHKJBGXX_K562_0min_celastrol10uM_rep2_GB_CAGATC_R1_minus.primary.bw");
file.bw.G8 <- c( "../k562/SRR182390x_plus.bw", "../k562/SRR182390x_minus.bw");
file.bw.GM <- c( "../GM12878/groseq_plus.bigWig", "../GM12878/groseq_minus.bigWig");
file.bw.GH <- c( "../HCT116/SRR1105736.7.plus.bw", "../HCT116/SRR1105736.7.minus.bw");
file.bw.CD4 <- c( "/fs/cbsudanko/storage/data/hg19/cd4/proseq/CD4-U_plus.bw", "/fs/cbsudanko/storage/data/hg19/cd4/proseq/CD4-U_minus.bw");
file.bw.MCF7 <- c( "/fs/cbsudanko/storage/data/hg19/mcf7/groseq/MCF7.unt.all_plus.bw", "/fs/cbsudanko/storage/data/hg19/mcf7/groseq/MCF7.unt.all_minus.bw");
file.bw.HELA <- c( "/fs/cbsudanko/storage/data/hg19/hela/groseq/groseq_plus.bigWig", "/fs/cbsudanko/storage/data/hg19/hela/groseq/groseq_minus.bigWig");

file.dREG.G1 <- "/workdir/zw355/proj/prj10-dreg/new-rf-201803/G1/G1.dREG.peak.score.bed.gz";
file.dREG.G2 <- "/workdir/zw355/proj/prj10-dreg/new-rf-201803/G2/G2.dREG.peak.score.bed.gz";
file.dREG.G3 <- "/workdir/zw355/proj/prj10-dreg/new-rf-201803/G3/G3.dREG.peak.score.bed.gz";
file.dREG.G4 <- "/workdir/zw355/proj/prj10-dreg/new-rf-201803/G4/G4.dREG.peak.score.bed.gz";
file.dREG.G5 <- "/workdir/zw355/proj/prj10-dreg/new-rf-201803/G5/G5.dREG.peak.score.bed.gz";
file.dREG.G6 <- "/workdir/zw355/proj/prj10-dreg/new-rf-201803/G6/G6.dREG.peak.score.bed.gz";
file.dREG.G7 <- "/workdir/zw355/proj/prj10-dreg/new-rf-201803/G7/G7.dREG.peak.score.bed.gz";
file.dREG.G8 <- "/workdir/zw355/proj/prj10-dreg/new-rf-201803/G8/G8.dREG.peak.score.bed.gz";
file.dREG.GM <- "/workdir/zw355/proj/prj10-dreg/new-rf-201803/GM/GM.dREG.peak.score.bed.gz";
file.dREG.GH <- "/workdir/zw355/proj/prj10-dreg/new-rf-201803/GH/GH.dREG.peak.score.bed.gz";
file.dREG.CD4 <- "/workdir/zw355/proj/prj10-dreg/new-rf-201803/CD4/CD4-U.dREG.peak.score.bed.gz";
file.dREG.MCF7 <- "/workdir/zw355/proj/prj10-dreg/new-rf-201803/MCF7/MCF7.unt.dREG.peak.score.bed.gz";
file.dREG.HELA <- "/workdir/zw355/proj/prj10-dreg/new-rf-201803/HELA/HELA.groseq.dREG.peak.score.bed.gz";

file.grocap="../k562/hg19.k562.new_hmm2b.post2.bed"

file.H3K27ac.K562.peak  <- "/fs/cbsudanko/storage/data/hg19/k562/histones/wgEncodeBroadHistoneK562H3k27acStdAln.bed.gz"
file.H3K27ac.GM.peak  <- "/fs/cbsudanko/storage/data/hg19/gm12878/histones/wgEncodeBroadHistoneGm12878H3k27acStdPk.broadPeak.gz"
file.H3K27ac.HCT116.peak <- "../HCT116/hg19.H3K27ac.narrow.peak.ENCFF450AGJ.bed"
file.H3K27ac.CD4.peak <- "/fs/cbsudanko/storage/data/hg19/cd4/epiRoadmap_histone/H3K27ac_peaks.narrowPeak"
file.H3K27ac.MCF7.peak <- "/fs/cbsudanko/storage/data/hg19/mcf7/histones/wgEncodeSydhHistoneMcf7H3k27acUcdPk.narrowPeak.gz"
file.H3K27ac.HELA.peak <- "/fs/cbsudanko/storage/data/hg19/hela/histones/wgEncodeBroadHistoneHelas3H3k27acStdPk.broadPeak.gz"

file.DHS.K562.peak    <- "./merge.broad.peak.k562.bed";
file.DHS.GM12878.peak <- "./merge.narrowPeak.gm.bed"
file.DHS.HCT116.peak  <- "./merge.narrowPeak.hct116.bed"
file.DHS.CD4.peak     <- "./merge.narrowPeak.cd4.bed"
file.DHS.MCF7.peak    <- "./merge.narrowPeak.mcf7.bed"
file.DHS.HELA.peak    <- "./merge.narrowPeak.hela.bed"


file.DHS.peak <- c(rep(file.DHS.K562.peak, 8), file.DHS.GM12878.peak, file.DHS.HCT116.peak, file.DHS.CD4.peak, file.DHS.MCF7.peak, file.DHS.HELA.peak );
file.H3K27ac.peak <- c(rep(file.H3K27ac.K562.peak, 8), file.H3K27ac.GM.peak, file.H3K27ac.HCT116.peak, file.H3K27ac.CD4.peak, file.H3K27ac.MCF7.peak, file.H3K27ac.HELA.peak );

file.dREG <-c( file.dREG.G1, file.dREG.G2, file.dREG.G3, file.dREG.G4, file.dREG.G5, file.dREG.G6, file.dREG.G7, file.dREG.G8, file.dREG.GM, file.dREG.GH, file.dREG.CD4, file.dREG.MCF7, file.dREG.HELA );
file.bws <- list(file.bw.G1, file.bw.G2, file.bw.G3, file.bw.G4, file.bw.G5, file.bw.G6, file.bw.G7, file.bw.G8, file.bw.GM, file.bw.GH, file.bw.CD4, file.bw.MCF7, file.bw.HELA );

make_big_bed<-function( file.nonneg,dist=100 )
{
	system( paste("cat ../k562/wgEncodeOpenChromDnaseK562PkV2.narrowPeak ../k562/GSM646567_hg19_wgEncodeUwDgfK562Pk.narrowPeak.txt ../k562/GSM646567_hg19_wgEncodeUwDgfK562Pk.macs2.narrowPeak | awk -v OFS='\\t' '{print $1,$2,$3}' - > ", file.nonneg ), intern=TRUE );

	nonneg_bed <- unique(read.table( pipe(paste("cat ", file.nonneg, file.grocap, " | sort-bed ", file.nonneg," | bedtools merge -i - -d 100 ", sep=" " ) ) )[,c(1:3)]);
	nonneg_bed[,2] <- nonneg_bed[,2] - dist
	idx.mis <- which(nonneg_bed[,2]<0);
	if(length(idx.mis)>0) nonneg_bed[idx.mis,2] <- 0;
	nonneg_bed[,3] <- nonneg_bed[,3] + dist
	write.table( nonneg_bed, file=file.nonneg, quote=F, row.name=F, col.names=F, sep="\t" );
}

make_GM12878_bed<-function( file.nonneg,dist=100 )
{
	nonneg_bed <- read.table("/fs/cbsudanko/storage/data/hg19/gm12878/dnase/uw.merge.narrowPeak.bed");
	nonneg_bed[,2] <- nonneg_bed[,2] - dist
	idx.mis <- which(nonneg_bed[,2]<0);
	if(length(idx.mis)>0) nonneg_bed[idx.mis,2] <- 0;
	nonneg_bed[,3] <- nonneg_bed[,3] + dist
	write.table( nonneg_bed, file=file.nonneg, quote=F, row.name=F, col.names=F, sep="\t" );
}

make_HCT116_bed<-function( file.nonneg,dist=100 )
{
	nonneg_bed <- read.table("../HCT116/hg19.dNase.HCT116.broad.peak.bed");
	nonneg_bed[,2] <- nonneg_bed[,2] - dist
	idx.mis <- which(nonneg_bed[,2]<0);
	if(length(idx.mis)>0) nonneg_bed[idx.mis,2] <- 0;
	nonneg_bed[,3] <- nonneg_bed[,3] + dist
	write.table( nonneg_bed, file=file.nonneg, quote=F, row.name=F, col.names=F, sep="\t" );
}

make_CD4_bed<-function( file.nonneg,dist=100 )
{
	nonneg_bed <- read.table("/fs/cbsudanko/storage/data/hg19/cd4/dnase/duke.cd4.merge.bed");
	nonneg_bed[,2] <- nonneg_bed[,2] - dist
	idx.mis <- which(nonneg_bed[,2]<0);
	if(length(idx.mis)>0) nonneg_bed[idx.mis,2] <- 0;
	nonneg_bed[,3] <- nonneg_bed[,3] + dist
	write.table( nonneg_bed, file=file.nonneg, quote=F, row.name=F, col.names=F, sep="\t" );
}


make_MCF7_bed<-function( file.nonneg,dist=100 )
{
	system( paste("zcat /fs/cbsudanko/storage/data/hg19/mcf7/dnase/mcf7.peaks.bed.gz /fs/cbsudanko/storage/data/hg19/mcf7/dnase/wgEncodeOpenChromDnaseMcf7Pk.narrowPeak.gz /fs/cbsudanko/storage/data/hg19/mcf7/dgf/Mcf7Estctrl0h.peaks_peaks.narrowPeak.gz  /fs/cbsudanko/storage/data/hg19/mcf7/dgf/GSM1024784_hg19_wgEncodeUwDnaseMcf7Est100nm1hPkRep1.narrowPeak.gz | awk -v OFS='\\t' '{print $1,$2,$3}' -  > ", file.nonneg ), intern=TRUE );

	nonneg_bed <- unique(read.table( pipe(paste("cat ", file.nonneg, file.grocap, " | sort-bed ", file.nonneg," | bedtools merge -i - -d 100 ", sep=" " ) ) )[,c(1:3)]);
	nonneg_bed[,2] <- nonneg_bed[,2] - dist
	idx.mis <- which(nonneg_bed[,2]<0);
	if(length(idx.mis)>0) nonneg_bed[idx.mis,2] <- 0;
	nonneg_bed[,3] <- nonneg_bed[,3] + dist
	write.table( nonneg_bed, file=file.nonneg, quote=F, row.name=F, col.names=F, sep="\t" );
}

make_HELA_bed<-function( file.nonneg,dist=100 )
{
	nonneg_bed <- read.table("/fs/cbsudanko/storage/data/hg19/hela/dnase/uw.merge.narrowPeak.bed");
	nonneg_bed[,2] <- nonneg_bed[,2] - dist
	idx.mis <- which(nonneg_bed[,2]<0);
	if(length(idx.mis)>0) nonneg_bed[idx.mis,2] <- 0;
	nonneg_bed[,3] <- nonneg_bed[,3] + dist
	write.table( nonneg_bed, file=file.nonneg, quote=F, row.name=F, col.names=F, sep="\t" );
}


make_big_bed( file.DHS.K562.peak );
make_GM12878_bed( file.DHS.GM12878.peak );
make_HCT116_bed( file.DHS.HCT116.peak );
make_CD4_bed( file.DHS.CD4.peak );
make_MCF7_bed( file.DHS.MCF7.peak );
make_HELA_bed( file.DHS.HELA.peak );

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

n.overlap <- n.novel <- n.total <- c();
for(i in 1:length(file.dREG ) )
{
	dREG_bed0 <- read.table( pipe(paste("bedtools subtract -A -a ", file.dREG[i], "-b", file.black.bed)));
	n.total <- c(n.total, NROW(dREG_bed0));
	
	file.pure.dREG = tempfile();
	write.table(dREG_bed0, file=file.pure.dREG, quote=F, row.names=F, col.names=F, sep="\t");

library(tools);
    if (file_ext(file.H3K27ac.peak[i])=="gz")
	    tb <- unique(read.table(pipe(paste("zcat ", file.H3K27ac.peak[i], " | bedtools intersect -a" , file.pure.dREG, "-b - -loj")))[,c(1:5)])
	else
	    tb <- unique(read.table(pipe(paste("cat ", file.H3K27ac.peak[i], " | bedtools intersect -a" , file.pure.dREG, "-b - -loj")))[,c(1:5)]);
	
	idx.novel <- which( unlist(apply(tb, 1, function(x){return(as.character(x[5])==".")})) );
	tb.nov <- tb [ idx.novel, 1:4]

	file1= tempfile();
	write.table(tb.nov, file=file1, quote=F, row.names=F, col.names=F, sep="\t");

	tb <- unique(read.table(pipe(paste("cat ", file.DHS.peak[i], " | bedtools intersect -a" , file1, "-b - -loj")))[,c(1:5)]);
	idx.novel <- which( unlist(apply(tb, 1, function(x){return(as.character(x[5])==".")})) );
	tb.nov <- tb[ idx.novel, 1:4];

	unlink(file1);
	unlink(file.pure.dREG);

	n.novel <- c( n.novel, NROW(tb.nov) );
	n.overlap <- c(n.overlap, NROW(tb)-NROW(tb.nov) );
}

unlink(file.DHS.peak);

basesCovered <- basesCovered[-4];
coverage <- coverage[-4];
datasize <- datasize[-4];
dataname <- dataname[-4]
n.total <- n.total[-4]
n.novel <- n.novel[-4]
n.overlap <- n.overlap[-4]
n.novel.ratio <- n.novel/n.total

cK562T <- "#7fc97f"
cK562  <- "#f0027f"
cGM    <- "#beaed4"
cHCT   <- "#fdc086"
cCD4   <- "#666666"
cMCF7  <- "#386cb0"
cHELA  <- "#bf5b17"

cols <- c( rep(cK562T, 7), cGM, cHCT, cCD4, cMCF7, cHELA);
cols[6] <- cols[7] <- cK562;

#"#86db63", "#86db63", "#86db63", "#86db63", "#86db63", "#f92367", "#f92367", "#5070fb", "#CCAA80" );

show(coverage);
#log.coverage <- log(coverage)
coverage <- coverage/1000/1000

show(n.novel);

if(1)
{
	pdf("fig-3A.pdf")

	plot( coverage, n.novel, main="", xlab="Number of mapped reads (M)", ylab="#Novel elements", ylim=c(0,max(n.novel)), pch=17, col=cols, cex=2)
	r <- summary(lm( n.novel~ log(coverage)));
	r.b <- r$coefficients[1,1]
	r.k <- r$coefficients[2,1]
	r.x <- seq( 1, 500, 1 );
	lines( r.x, r.b+r.k*log(r.x), lty=22, lwd=1, col="black");

	legend.name <- c( "K562 training", "K562 holdout", "GM12878 holdout", "HCT116 holdout", "CD4 holdout", "MCF7 holdout", "HELA holdout");
	legend.cols  <- c( cK562T, cK562, cGM, cHCT, cCD4, cMCF7, cHELA );
	legend.pch  <- rep(17,7);
	legend("bottomright", legend.name, text.col=legend.cols, col=legend.cols, pch=legend.pch, cex=1, ncol=1, horiz=F);

	dev.off();
}

if(1)
{
	pdf("fig-3A2.pdf")

	n.novel.ratio[5] <- n.novel.ratio[5];
	plot( coverage, n.novel.ratio, main="", xlab="Number of mapped reads (M)", ylab="Fractions", ylim=c(0,max(n.novel.ratio)*1.05), pch=17, col=cols, cex=2)
	r <- summary(lm( n.novel.ratio~ log(coverage)));
	r.b <- r$coefficients[1,1]
	r.k <- r$coefficients[2,1]

	r.x <- seq( 1, 500, 1 );
	lines( r.x, r.b+r.k*log(r.x), lty=22, lwd=1, col="black");

	text(coverage[1], n.novel.ratio[1]+0.005, "G1", cex=0.8);
	text(coverage[2], n.novel.ratio[2]+0.005, "G2", cex=0.8);
	text(coverage[3], n.novel.ratio[3]+0.005, "G3", cex=0.8);
	text(coverage[4], n.novel.ratio[4]+0.005, "G5", cex=0.8);
	text(coverage[5], n.novel.ratio[5]-0.004, "G6", cex=0.8);
	text(coverage[6]*1.25, n.novel.ratio[6]+0.003, "G7", cex=0.8);
	text(coverage[7], n.novel.ratio[7]+0.005, "G8", cex=0.8);
	text(coverage[8], n.novel.ratio[8]+0.005, "GM", cex=0.8);
	text(coverage[9], n.novel.ratio[9]+0.005, "GH", cex=0.8);
	text(coverage[10], n.novel.ratio[10]+0.005, "CD4", cex=0.8);
	text(coverage[11], n.novel.ratio[11]+0.005, "MCF7", cex=0.8);
	text(coverage[12], n.novel.ratio[12]+0.005, "H1", cex=0.8);


	legend.name <- c( "K562 training", "K562 holdout", "GM12878 holdout", "HCT116 holdout", "CD4 holdout", "MCF7 holdout", "HELA holdout");
	legend.cols  <- c( cK562T, cK562, cGM, cHCT, cCD4, cMCF7, cHELA );
	legend.pch  <- rep(17,7);
	legend("bottomright", legend.name, text.col=legend.cols, col=legend.cols, pch=legend.pch, cex=1, ncol=1, horiz=F);

	dev.off();
}