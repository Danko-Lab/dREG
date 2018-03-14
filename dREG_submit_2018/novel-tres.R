library(bigWig);


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
file.K562.DHS.peak  <- "../k562/K562.bed";

file.GM.H3K27ac.peak  <- "/fs/cbsudanko/storage/data/hg19/gm12878/histones/wgEncodeBroadHistoneGm12878H3k27acStdPk.broadPeak.gz"

file.dREG.G1 <- "/workdir/zw355/proj/prj10-dreg/new-model-201712/G1/G1.bedgraph.dREG.peak.score.bed.gz";
file.dREG.G2 <- "/workdir/zw355/proj/prj10-dreg/new-model-201712/G2/G2.bedgraph.dREG.peak.score.bed.gz";
file.dREG.G3 <- "/workdir/zw355/proj/prj10-dreg/new-model-201712/G3/G3.bedgraph.dREG.peak.score.bed.gz";
file.dREG.G4 <- "/workdir/zw355/proj/prj10-dreg/new-model-201712/G4/G4.bedgraph.dREG.peak.score.bed.gz";
file.dREG.G5 <- "/workdir/zw355/proj/prj10-dreg/new-model-201712/G5/G5.bedgraph.dREG.peak.score.bed.gz";
file.dREG.G6 <- "/workdir/zw355/proj/prj10-dreg/new-model-201712/G6/G6.bedgraph.dREG.peak.score.bed.gz";
file.dREG.G7 <- "/workdir/zw355/proj/prj10-dreg/new-model-201712/G7/G7.bedgraph.dREG.peak.score.bed.gz";
file.dREG.GM <- "/workdir/zw355/proj/prj10-dreg/new-model-201712/GM/GM.bedgraph.dREG.peak.score.bed.gz";
file.dREG <-c( file.dREG.G1, file.dREG.G2, file.dREG.G3, file.dREG.G4, file.dREG.G5, file.dREG.G6, file.dREG.G7, file.dREG.GM);


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


dist_X_range <- function( pred.bed, file.X )
{
	options("scipen"=100, "digits"=4);

	tmp.bed <- tempfile(fileext=".bed" )
	write.table( pred.bed[,c(1:3)], file=tmp.bed, row.names=F, col.names=F, quote=F, sep="\t");
	system( paste("sort-bed ",  tmp.bed, " > ", tmp.bed, ".sorted", sep="") );

	tb.close <-  read.table( file = pipe(paste("sort-bed ",  file.X, " | bedtools closest -d -a ", tmp.bed, ".sorted -b - -t first", sep="") ) );
	return(tb.close);
}


i<-1;
tb <- unique(read.table(pipe(paste("zcat ", file.K562.H3K27ac.peak , " | bedtools intersect -a" , file.dREG[i], "-b - -loj")))[,c(1:5)]);
idx.novel <- which( unlist(apply(tb, 1, function(x){return(as.character(x[5])==".")})) );
tb.nov <- tb [ idx.novel, 1:4]
		
file1= tempfile();
write.table(tb.nov, file=file1, quote=F, row.names=F, col.names=F, sep="\t");

tb <- unique(read.table(pipe(paste("cat ", file.K562.DHS.peak, " | bedtools intersect -a" , file1, "-b - -loj")))[,c(1:5)]);
idx.novel <- which( unlist(apply(tb, 1, function(x){return(as.character(x[5])==".")})) );
tb.nov <- tb[ idx.novel, 1:4];
		
unlink(file1);

tb.center <- dist_X_center(tb.nov, file.K562.DHS.peak);
tb.range <- dist_X_range(tb.nov, file.K562.DHS.peak);


