library(bigWig);

# H3k27ac:
# K562: /fs/cbsudanko/storage/data/hg19/k562/histones/wgEncodeBroadHistoneK562H3k27acStdAln.bed.gz
TPR.H3k27ac <- c(G1=0.786, G2=0.526, G3=0.733, G4=0.635, G5=0.687, G6=0.582, G7=0.584, G8=0.654, GM=0.718, GH=0.849, CD4=0.814, MCF7=0.842, HELA=0.450 )
#../k562/hg19.k562.grocap.pair.bed
TPR.grocap <-  c(G1=0.949, G2=0.782, G3=0.923, G4=0.840, G5=0.894, G6=0.811, G7=0.813, G8=0.819, GM=0.961, NA, NA, NA, NA)

dataname <- c( paste("K562-G", 1:8, sep=""),         "GM12878", "HCT116", "CD4",     "MCF7",    "HELA" )

#http://colorbrewer2.org/#type=qualitative&scheme=Accent&n=7
cK562T <- "#7fc97f"
cK562  <- "#f0027f"
cGM    <- "#beaed4"
cHCT   <- "#fdc086"
cCD4   <- "#666666"
cMCF7  <- "#386cb0"
cHELA  <- "#bf5b17"
cols <-     c( rep(cK562T,6), cK562, cK562, cGM, cHCT, cCD4, cMCF7, cHELA );


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

file.bws <- list(file.bw.G1, file.bw.G2, file.bw.G3, file.bw.G4, file.bw.G5, file.bw.G6, file.bw.G7, file.bw.G8, file.bw.GM, file.bw.GH, file.bw.CD4, file.bw.MCF7, file.bw.HELA);
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

TPR.H3k27ac <- TPR.H3k27ac[-4];
TPR.grocap <- TPR.grocap[-4];
basesCovered <- basesCovered[-4];
coverage <- coverage[-4];
datasize <- datasize[-4];
dataname <- dataname[-4];
cols <- cols[-4];

show(coverage);

coverage <-coverage/1000/1000

pdf("fig-1C.pdf")

plot( coverage, TPR.grocap, main="", xlab="Number of mapped reads (M)", ylab="Sensitivity", pch=19, col=cols, ylim=c(0.45, 1), cex=2)
r <- summary(lm( TPR.grocap~log(coverage)));
r.b <- r$coefficients[1,1]
r.k <- r$coefficients[2,1]
r.x <- seq( 0, 500, 1);
lines( r.x, r.b+r.k*log(r.x), lty=22, lwd=1, col="black");

points(coverage, TPR.H3k27ac, pch=17, col=cols, cex=2)
r <- summary(lm( TPR.H3k27ac~log(coverage)));
r.b <- r$coefficients[1,1]
r.k <- r$coefficients[2,1]
#r.x <- c( seq( min(log.coverage), max(log.coverage), log(1000000*0.1) ),  max(log.coverage))
r.x <- seq( 0, 500, 1);
lines( r.x, r.b+r.k*log(r.x), lty=22, lwd=1, col="black");

#axis(1, at=seq(16.5, 20, 0.7), labels=round( exp(seq(16.5,20, 0.7))/1000/1000, 1), cex=2);

#legend.name <- c("Training/GRO-cap", "Training/H3k27ac", "K562 holdout/GRO-cap", "K562 holdout/H3k27ac", "GM12878 holdout/GRO-cap", "GM12878 holdout/H3k27ac");
#legend.cols  <- c("#86db63", "#86db63", "#f92367", "#f92367", "#5070fb", "#5070fb")
#legend.pch  <- c(19, 17, 19, 17, 19, 17)
#legend("bottomright", legend.name[c(1,3,5,2,4,6)], text.col=legend.cols[c(1,3,5,2,4,6)], col=legend.cols[c(1,3,5,2,4,6)], pch=legend.pch[c(1,3,5,2,4,6)] , cex=1, ncol=2, horiz=F);

legend.name <- c("GRO-cap Pairs", "Training", "K562 holdout", "GM12878 holdout", "HCT116 holdout", "", "", "",
                 "H3k27ac/DHS", "Training", "K562 holdout", "GM12878 holdout", "HCT116 holdout", "CD4 holdout", "MCF7 holdout", "HELA holdout");

legend.cols  <- c("black", cK562T, cK562, cGM, "white", "white", "white", "white", 
                  "black", cK562T, cK562, cGM, cHCT,    cCD4,    cMCF7,    cHELA )
legend.pch  <- c( -1, 19, 19, 19, NA, NA,NA,NA,
                  -1, 17, 17, 17, 17, 17, 17, 17)
legend("bottomright", legend.name, text.col=legend.cols, col=legend.cols, pch=legend.pch, cex=1, ncol=2, horiz=F);

dev.off();

