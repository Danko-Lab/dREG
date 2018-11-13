library(bigWig);

#/fs/cbsudanko/storage/data/hg19/k562/histones/wgEncodeBroadHistoneK562H3k27acStdAln.bed.gz
TPR.H3k27ac <- c(0.8027,0.5486,0.7470,0.6707,0.7035,0.5959,0.6006,0.8346)
#../k562/hg19.k562.grocap.pair.bed
TPR.grocap <- c(0.963, 0.820, 0.940, 0.837, 0.916, 0.840, 0.843, 0.970)

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

file.TF.chipseq <- "/fs/cbsudanko/storage/data/hg19/all/ENCODE_tf_peak_calls/wgEncodeRegTfbsClusteredWithCellsV3.bed.gz"

#file.K562.DHS.peak  <- "../k562/K562.bed";
file.K562.DHS.peak  <- "./k562.merged.broad.peak.bed";
file.grocap="../k562/hg19.k562.new_hmm2b.post2.bed"
file.K562.H3K27ac.peak  <- "/fs/cbsudanko/storage/data/hg19/k562/histones/wgEncodeBroadHistoneK562H3k27acStdAln.bed.gz"
file.K562.H3kme1.peak   <- "/fs/cbsudanko/storage/data/hg19/k562/histones/wgEncodeBroadHistoneK562H3k4me1StdAln.bed.gz"


file.GM12878.H3K27ac.peak  <- "/fs/cbsudanko/storage/data/hg19/gm12878/histones/wgEncodeBroadHistoneGm12878H3k27acStdPk.broadPeak.gz"
file.GM12878.DHS.peak <- "/fs/cbsudanko/storage/data/hg19/gm12878/dnase/uw.merge.narrowPeak.bed"
file.GM12878.H3kme1.peak <- "/fs/cbsudanko/storage/data/hg19/gm12878/histones/wgEncodeBroadHistoneGm12878H3k4me1StdPk.broadPeak.gz"

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

options("scipen"=100, "digits"=4);

n.overlap <- n.novel <- n.less500 <- n.less1k <- c();
for(i in 1:length(file.dREG ) )
{
    file.H3K27ac.peak <- file.K562.H3K27ac.peak;
    file.DHS.peak <- file.K562.DHS.peak;
    file.H3kme1.peak<-  file.K562.H3kme1.peak;
	if(i==8)
	{
	    file.H3K27ac.peak <- file.GM12878.H3K27ac.peak;
	    file.DHS.peak <- file.GM12878.DHS.peak;
	    file.H3kme1.peak<-  file.GM12878.H3kme1.peak;
	}

	tb <- unique(read.table(pipe(paste("zcat ", file.H3K27ac.peak , " | bedtools intersect -a" , file.dREG[i], "-b - -loj")))[,c(1:5)]);
	idx.novel <- which( unlist(apply(tb, 1, function(x){return(as.character(x[5])==".")})) );
	tb.nov <- tb [ idx.novel, 1:4]
		
	file1= tempfile();
	write.table(tb.nov, file=file1, quote=F, row.names=F, col.names=F, sep="\t");

	tb <- unique(read.table(pipe(paste("cat ", file.DHS.peak, " | bedtools intersect -a" , file1, "-b - -loj")))[,c(1:5)]);
	idx.novel <- which( unlist(apply(tb, 1, function(x){return(as.character(x[5])==".")})) );
	tb.nov <- tb[ idx.novel, 1:4];

	write.table(tb.nov, file=file1, quote=F, row.names=F, col.names=F, sep="\t");
	tb.close <-  read.table( file = pipe(paste("sort-bed ", file1, " | bedtools closest -d -a  - -b ",file.DHS.peak, " -t first", sep="") ) );
	
	n.novel <- c(n.novel, NROW(tb.nov));
	n.less500 <- c(n.less500,  sum(tb.close$V8<500) );
	n.less1k <- c(n.less1k,  sum(tb.close$V8<1000) - sum(tb.close$V8<500) );

	unlink(file1);

}

cols <- c("#86db63", "#86db63", "#86db63", "#86db63", "#86db63", "#f92367", "#5070fb" );

show(n.novel);

dat <- rbind( n.less500[-4], n.less1k[-4], (n.novel - n.less500 - n.less1k)[-4] );
mat <- t(as.matrix(t(dat)/colSums(dat)))
colnames(mat) <- c("G1", "G2", "G3", "G5", "G6", "G7", "GM");

pdf("fig-S9.pdf")
barplot(mat, xlab="Novel elements", ylab="<500(drak), <1000(medium), >=1000(light)" );
dev.off();

