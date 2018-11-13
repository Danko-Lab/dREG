library(bigWig);

dataname <- c(paste("K562-G", 1:8, sep=""), "GM12878", "HCT116")

file.black.bed <- "/fs/cbsudanko/storage/data/hg19/all/encode_blacklist/hg19.blacklist.bed";

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

#file.K562.DHS.peak  <- "../k562/K562.bed";
file.K562.DHS.peak  <- "./k562.merged.broad.peak.bed";
file.grocap="../k562/hg19.k562.new_hmm2b.post2.bed"
file.K562.H3K27ac.peak  <- "/fs/cbsudanko/storage/data/hg19/k562/histones/wgEncodeBroadHistoneK562H3k27acStdAln.bed.gz"
file.K562.H3kme1.peak   <- "/fs/cbsudanko/storage/data/hg19/k562/histones/wgEncodeBroadHistoneK562H3k4me1StdAln.bed.gz"

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


extract_TF <- function(celline, file.TF.chipseq.cell)
{
	file.TF.chipseq <- "/fs/cbsudanko/storage/data/hg19/all/ENCODE_tf_peak_calls/wgEncodeRegTfbsClusteredWithCellsV3.bed.gz"
	system(paste( "zcat", file.TF.chipseq, " | grep ", celline, " | cut -f1,2,3 > ", file.TF.chipseq.cell ) );
}

file.TF.chipseq.K562 <- "wgEncodeRegTfbsClusteredWithCellsV3.K562.bed.gz"
extract_TF( "K562", file.TF.chipseq.K562 );

file.TF.chipseq.GM12878 <- "wgEncodeRegTfbsClusteredWithCellsV3.GM12878.bed.gz"
extract_TF( "GM12878", file.TF.chipseq.GM12878 );

file.TF.chipseq.HCT116 <- "wgEncodeRegTfbsClusteredWithCellsV3.HCT116.bed.gz"
extract_TF( "HCT-116", file.TF.chipseq.HCT116 );

file.GM12878.DHS.peak <- "/fs/cbsudanko/storage/data/hg19/gm12878/dnase/uw.merge.narrowPeak.bed"
file.GM12878.H3K27ac.peak  <- "/fs/cbsudanko/storage/data/hg19/gm12878/histones/wgEncodeBroadHistoneGm12878H3k27acStdPk.broadPeak.gz"
file.GM12878.H3kme1.peak <- "/fs/cbsudanko/storage/data/hg19/gm12878/histones/wgEncodeBroadHistoneGm12878H3k4me1StdPk.broadPeak.gz"

file.HCT116.DHS.peak <- "../HCT116/hg19.dNase.HCT116.broad.peak.bed"
file.HCT116.H3K27ac.peak  <- "../HCT116/hg19.H3K27ac.narrow.peak.ENCFF450AGJ.bed.gz"
file.HCT116.H3kme1.peak <- "../HCT116/hg19.H3K4me1.narrow.peak.ENCFF884XUC.bed.gz"

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
file.dREG <-c( file.dREG.G1, file.dREG.G2, file.dREG.G3, file.dREG.G4, file.dREG.G5, file.dREG.G6, file.dREG.G7, file.dREG.G8, file.dREG.GM, file.dREG.GH );

file.bws <- list(file.bw.G1, file.bw.G2, file.bw.G3, file.bw.G4, file.bw.G5, file.bw.G6, file.bw.G7, file.bw.G8, file.bw.GM, file.bw.GH );
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

n.overlap <- n.novel <- n.removeTFBS <- n.removeH3k4me1 <- n.DHS <- n.H3K27ac <- n.total <- nall.removeTFBS <- nall.removeH3k4me1 <-c();
for(i in 1:length(file.dREG ) )
{
    file.H3K27ac.peak <- file.K562.H3K27ac.peak;
    file.DHS.peak <- file.K562.DHS.peak;
    file.H3kme1.peak<-  file.K562.H3kme1.peak;
    file.TF.chipseq <- file.TF.chipseq.K562
	if(i==9)
	{
	    file.H3K27ac.peak <- file.GM12878.H3K27ac.peak;
	    file.DHS.peak <- file.GM12878.DHS.peak;
	    file.H3kme1.peak<-  file.GM12878.H3kme1.peak;
	    file.TF.chipseq <- file.TF.chipseq.GM12878
	}
    else
	if(i==10)
	{
	    file.H3K27ac.peak <- file.HCT116.H3K27ac.peak;
	    file.DHS.peak <- file.HCT116.DHS.peak;
	    file.H3kme1.peak<-  file.HCT116.H3kme1.peak;
	    file.TF.chipseq <- file.TF.chipseq.HCT116
	}

    #remove black list
	dREG_bed0 <- read.table( pipe(paste("bedtools subtract -A -a ", file.dREG[i], "-b", file.black.bed)));
	n.total <- c(n.total, NROW(dREG_bed0));
	
	file.pure.dREG = tempfile();
	write.table(dREG_bed0, file=file.pure.dREG, quote=F, row.names=F, col.names=F, sep="\t");

    #intersect with DHS
	tb <- unique(read.table(pipe(paste("cat ", file.DHS.peak , " | bedtools intersect -a" , file.pure.dREG, "-b - -loj")))[,c(1:5)]);
	idx.novel <- which( unlist(apply(tb, 1, function(x){return(as.character(x[5])==".")})) );
	n.DHS <- c(n.DHS, NROW(tb) - NROW(idx.novel));
	tb.nov <- tb [ idx.novel, 1:4]

	file1= tempfile();
	write.table(tb.nov, file=file1, quote=F, row.names=F, col.names=F, sep="\t");

    #intersect with H3k27ac
	tb <- unique(read.table(pipe(paste("zcat ", file.H3K27ac.peak, " | bedtools intersect -a" , file1, "-b - -loj")))[,c(1:5)]);
	idx.novel <- which( unlist(apply(tb, 1, function(x){return(as.character(x[5])==".")})) );
	n.H3K27ac <- c(n.H3K27ac, NROW(tb) - NROW(idx.novel));
	tb.nov <- tb[ idx.novel, 1:4];

	n.novel <- c( n.novel, NROW(tb.nov) );
	n.overlap <- c(n.overlap, NROW(tb)-NROW(tb.nov) );

    #intersect with TF chipseq
	write.table(tb.nov, file=file1, quote=F, row.names=F, col.names=F, sep="\t");
	tb <- unique(read.table(pipe(paste("cat ", file.TF.chipseq, " | bedtools intersect -a" , file1, "-b - -loj")))[,c(1:5)])
	idx.novel <- which( unlist(apply(tb, 1, function(x){return(as.character(x[5])==".")})) );
	tb.nov2 <- tb[ idx.novel, 1:4];

	n.removeTFBS <- c( n.removeTFBS, NROW(tb.nov) - NROW(tb.nov2) );

	tb.nov <- tb.nov2;
	write.table(tb.nov, file=file1, quote=F, row.names=F, col.names=F, sep="\t");

    #intersect with H3k4me1
	tb <- unique(read.table(pipe(paste("cat ", file.H3kme1.peak, " | bedtools intersect -a" , file1, "-b - -loj")))[,c(1:5)]);
	idx.novel <- which( unlist(apply(tb, 1, function(x){return(as.character(x[5])==".")})) );
	tb.nov2 <- tb[ idx.novel, 1:4];

	n.removeH3k4me1<- c( n.removeH3k4me1, NROW(tb.nov) - NROW(tb.nov2) );

    # calcuate same count based on all elemenents
    #intersect with TF chipseq
    tb.nov <- dREG_bed0;
	write.table(tb.nov, file=file1, quote=F, row.names=F, col.names=F, sep="\t");
	tb <- unique(read.table(pipe(paste("cat ", file.TF.chipseq, " | bedtools intersect -a" , file1, "-b - -loj")))[,c(1:5)])
	idx.novel <- which( unlist(apply(tb, 1, function(x){return(as.character(x[5])==".")})) );
	tb.nov2 <- tb[ idx.novel, 1:4];

	nall.removeTFBS <- c( nall.removeTFBS, NROW(tb.nov) - NROW(tb.nov2) );

    #intersect with H3kme1
	tb.nov <- tb.nov2;
	write.table(tb.nov, file=file1, quote=F, row.names=F, col.names=F, sep="\t");
	tb <- unique(read.table(pipe(paste("cat ", file.H3kme1.peak, " | bedtools intersect -a" , file1, "-b - -loj")))[,c(1:5)]);
	idx.novel <- which( unlist(apply(tb, 1, function(x){return(as.character(x[5])==".")})) );
	tb.nov2 <- tb[ idx.novel, 1:4];

	nall.removeH3k4me1<- c( nall.removeH3k4me1, NROW(tb.nov) - NROW(tb.nov2) );

	unlink(file1);
	
	unlink(file.pure.dREG);

}

basesCovered <- basesCovered[-4];
coverage <- coverage[-4];
datasize <- datasize[-4];
dataname <- dataname[-4]

cols <- c("#86db63", "#86db63", "#86db63", "#86db63", "#86db63", "#f92367", "#5070fb" );

show(coverage);
#log.coverage <- log(coverage)
coverage <- coverage/1000/1000

show(n.novel);

mat0 <- rbind( n.removeTFBS[-4], n.removeH3k4me1[-4], (n.novel - n.removeTFBS - n.removeH3k4me1)[-4] );
colnames(mat0) <- c("G1", "G2", "G3", "G5", "G6", "G7", "G8","GM", "GH");
pdf("fig-S8.pdf")
barplot(mat0, xlab="Novel elements", ylab="#Novel elements", legend = c("Overlapping transcription sites", "Overlapping H3k4me1", "Others"), args.legend = list(x = "top")  );
dev.off();


pdf("fig-S8a.pdf")
mat1 <- rbind( nall.removeTFBS[-4], nall.removeH3k4me1[-4], (n.total - nall.removeTFBS - nall.removeH3k4me1)[-4] );
colnames(mat1) <- c("G1", "G2", "G3", "G5", "G6", "G7", "G8","GM", "GH");
barplot(mat1, xlab="All elements", ylab="#All elements", legend = c("Overlapping transcription sites", "Overlapping H3k4me1", "Others"), args.legend = list(x = "top") );
dev.off();


pdf("fig-S8b.pdf")
mat2 <- rbind( n.DHS[-4], n.H3K27ac[-4], n.novel[-4] );
colnames(mat2) <- c("G1", "G2", "G3", "G5", "G6", "G7", "G8","GM", "GH");
barplot(mat2, xlab="All elements", ylab="#All elements", legend = c("DHS overlapped", "H3K27ac overlapped", "Novel elements"), args.legend = list(x = "top") );
dev.off();

mat <- cbind(t(mat0), t(mat1), t(mat2))
colnames(mat) <- c("O_TS_novel", "O_H3k4me1_novel", "O_others_novel", "O_TS_all", "O_H3k4me1_all", "O_others_all","O_DHS_all", "O_H3k27ac_all", "novel_all" );

options(width=250);
show(mat)

unlink(file.TF.chipseq.K562)
unlink(file.TF.chipseq.GM12878)
unlink(file.TF.chipseq.HCT116)
unlink(file.K562.DHS.peak);

