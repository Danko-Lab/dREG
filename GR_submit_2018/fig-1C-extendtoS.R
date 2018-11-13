library(bigWig);

# H3k27ac:
# K562: /fs/cbsudanko/storage/data/hg19/k562/histones/wgEncodeBroadHistoneK562H3k27acStdAln.bed.gz
TPR.H3k27ac <- c(G1=0.786, G2=0.526, G3=0.733, G4=0.635, G5=0.687, G6=0.582, G7=0.584, G8=0.654, GM=0.718, GH=0.849, CD4=0.814, MCF7=0.842, HELA=0.450,  Churchman=0.639, Proudfoot=0.660)
EXP.Type <- c(G1="P", G2="G", G3="P", G4="P", G5="P", G6="P", G7="P", G8="G", GM="G", GH="G", CD4="P", MCF7="G", HELA="G", Churchman="N", Proudfoot="N")

dataname <- c( paste("K562-G", 1:8, sep=""),         "GM12878", "HCT116", "CD4",     "MCF7",    "HELA", "HELA", "HELA" )

#http://colorbrewer2.org/#type=qualitative&scheme=Accent&n=7
cPROSeq  <- "#7fc97f"
cGROSeq  <- "#f0027f"
ccNETSeq  <- "#beaed4"
cpNETSeq   <- "#386cb0"
#cpNETSeq    <- "#fdc086"
#cCD4    <- "#666666"
#cHELA   <- "#bf5b17"
cols <-     c( G1=cPROSeq, G2=cGROSeq, G3=cPROSeq, G4=cPROSeq, G5=cPROSeq, G6=cPROSeq, G7=cPROSeq, G8=cGROSeq, GM=cGROSeq, GH=cGROSeq, CD4=cPROSeq, MCF7=cGROSeq, HELA=cGROSeq,  Churchman=ccNETSeq, Proudfoot=cpNETSeq );


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

file.bw.Churchman <- c( "../mNET-seq/HeLaS3_Rep12.plus.bw", "../mNET-seq/HeLaS3_Rep12.minus.bw");
file.bw.Proudfoot <- c( "../mNET-seq/GSM1474225_ANET_8WG16_rep1_F.bw", "../mNET-seq/GSM1474225_ANET_8WG16_rep1_R.bw");

file.bws <- list(file.bw.G1, file.bw.G2, file.bw.G3, file.bw.G4, file.bw.G5, file.bw.G6, file.bw.G7, file.bw.G8, file.bw.GM, file.bw.GH, file.bw.CD4, file.bw.MCF7, file.bw.HELA, file.bw.Churchman, file.bw.Proudfoot);
coverage <- basesCovered <- datasize <- c();
for(i in 1:length(file.bws))
{
	bw1 <- load.bigWig(file.bws[[i]][1]);
	bw2 <- load.bigWig(file.bws[[i]][2]);
	show(str(bw1));
	coverage <-c(coverage, abs(bw1$mean*bw1$basesCovered) + abs(bw2$mean*bw2$basesCovered));
	basesCovered<-c(basesCovered, bw1$basesCovered+bw2$basesCovered );
	datasize<-c(datasize, bw1$primaryDataSize+bw2$primaryDataSize );
	unload.bigWig(bw1);
	unload.bigWig(bw2);
}

EXP.Type  <- EXP.Type[-4];
TPR.H3k27ac <- TPR.H3k27ac[-4];
basesCovered <- basesCovered[-4];
coverage <- coverage[-4];
datasize <- datasize[-4];
dataname <- dataname[-4];
cols <- cols[-4];

show(coverage);

coverage <-coverage/1000/1000

pdf("fig-1C-extendtoS.pdf")

plot( coverage, TPR.H3k27ac, main="", xlab="Number of mapped reads (M)", ylab="Fraction of Sites Recovered", pch=19, col=cols, ylim=c(0.45, 1), cex=2)
r <- summary(lm( TPR.H3k27ac~log(coverage)));
r.b <- r$coefficients[1,1]
r.k <- r$coefficients[2,1]
r.x <- seq( 0, 500, 1);
lines( r.x, r.b+r.k*log(r.x), lty=22, lwd=1, col="black");


legend.name <- c("PRO-seq", "GRO-seq", "mNET-seq(Churchman)", "mNET-seq(Proudfoot )");
legend.cols  <- c(cPROSeq, cGROSeq, ccNETSeq, cpNETSeq)
legend.pch  <- c( 19, 19, 19, 19 )
legend("bottomright", legend.name, text.col=legend.cols, col=legend.cols, pch=legend.pch, cex=1);

dev.off();

show(data.frame( EXP.Type, TPR.H3k27ac, coverage));