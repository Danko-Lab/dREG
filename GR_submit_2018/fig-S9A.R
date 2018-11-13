library(bigWig);

dataname <- c(paste("K562-G", 1:8, sep=""), "GM12878", "HCT116")

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


file.black.bed <- "/fs/cbsudanko/storage/data/hg19/all/encode_blacklist/hg19.blacklist.bed";
file.grocap="../k562/hg19.k562.new_hmm2b.post2.bed"
file.K562.H3K27ac.peak  <- "/fs/cbsudanko/storage/data/hg19/k562/histones/wgEncodeBroadHistoneK562H3k27acStdAln.bed.gz"
file.K562.H3kme1.peak   <- "/fs/cbsudanko/storage/data/hg19/k562/histones/wgEncodeBroadHistoneK562H3k4me1StdAln.bed.gz"

file.GM12878.DHS.peak <- "/fs/cbsudanko/storage/data/hg19/gm12878/dnase/uw.merge.narrowPeak.bed"
file.GM12878.H3K27ac.peak  <- "/fs/cbsudanko/storage/data/hg19/gm12878/histones/wgEncodeBroadHistoneGm12878H3k27acStdPk.broadPeak.gz"
file.GM12878.H3kme1.peak <- "/fs/cbsudanko/storage/data/hg19/gm12878/histones/wgEncodeBroadHistoneGm12878H3k4me1StdPk.broadPeak.gz"

file.HCT116.DHS.peak <- "../HCT116/hg19.dNase.HCT116.broad.peak.bed"
file.HCT116.H3K27ac.peak  <- "../HCT116/hg19.H3K27ac.narrow.peak.ENCFF450AGJ.bed.gz"
file.HCT116.H3kme1.peak <- "../HCT116/hg19.H3K4me1.narrow.peak.ENCFF884XUC.bed.gz"

#n.nodREG_Cap <- n.dREG_NoCap <- n.dREG_Cap <- c();
#for(i in 1:length(file.dREG ) )
dREG_bed0 <- read.table( pipe(paste("bedtools subtract -A -a ", file.dREG.G1, "-b", file.black.bed)));
n.total <- NROW(dREG_bed0);
	
file.pure.dREG = tempfile();
write.table(dREG_bed0, file=file.pure.dREG, quote=F, row.names=F, col.names=F, sep="\t");

tb <- unique(read.table(pipe(paste("bedtools intersect -a" , file.grocap, "-b", file.pure.dREG, "-loj"))))
idx.nomatch <- which(as.character(tb[,7])==".");
n.cap_nodREG <- NROW(idx.nomatch);
n.cap_dREG <- NROW(tb) - NROW(idx.nomatch);

tb <- unique(read.table(pipe(paste("bedtools intersect -a" , file.pure.dREG, "-b ",  file.grocap, "-loj"))))
idx.nomatch <- which(as.character(tb[,5])==".");
n.dREG_Cap <- NROW(tb) - NROW(idx.nomatch);
n.dREG_NoCap <- NROW(idx.nomatch);

unlink(file.pure.dREG);

show(n.dREG_Cap==n.cap_dREG);

if(0)
{
  library(venn);
  pdf("fig-S9A.pdf");
  x <- array(0, dim=c(n.dREG_NoCap + n.cap_nodREG + n.dREG_Cap, 2));
  x[1:(n.dREG_NoCap + n.dREG_Cap),1] <- 1;
  x[,2] <- 1;
  x[1:n.dREG_NoCap,2] <- 0;
  colnames(x) <- c("dREG","GRO-cap");

  venn(as.data.frame(x),snames=c("dREG", "GRO-cap"), zcolor="red, green");
  dev.off();
}

if(1)
{
  tb.dREG.noCap <- tb[idx.nomatch,1:6]
  tb.dREG.Cap <- tb[-idx.nomatch,1:6]

  tb.dREG.noCap.100 <- tb.dREG.noCap;
  tb.dREG.noCap.100[,3] <- tb.dREG.noCap.100[,2] +100
  tb.dREG.Cap.100 <- tb.dREG.Cap;
  tb.dREG.Cap.100[,3] <- tb.dREG.Cap.100[,2] +100

  library(bigWig);
  G1_pl <- load.bigWig( file.bw.G1[1] );
  G1_mn <- load.bigWig( file.bw.G1[2] );

  reads.noCap <- bed.region.bpQuery.bigWig(G1_pl, tb.dREG.noCap.100[,c(1:3)], abs.value = TRUE) + bed.region.bpQuery.bigWig(G1_mn, tb.dREG.noCap.100[,c(1:3)], abs.value = TRUE);
  reads.Cap   <- bed.region.bpQuery.bigWig(G1_pl, tb.dREG.Cap.100[,c(1:3)], abs.value = TRUE)   + bed.region.bpQuery.bigWig(G1_mn, tb.dREG.Cap.100[,c(1:3)], abs.value = TRUE);
               
  pdf("fig-S9-x.pdf");
  hist(log(reads.noCap), breaks=100, xlim=c(0, 10));
  hist(log(reads.Cap), breaks=100, xlim=c(0, 10));
  dev.off()

}




