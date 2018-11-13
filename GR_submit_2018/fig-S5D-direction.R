
file.DREG.G1 <- "../new-rf-201803/G1/G1.dREG.peak.full.bed.gz";
file.DREG.G2 <- "../new-rf-201803/G2/G2.dREG.peak.full.bed.gz";
file.DREG.G3 <- "../new-rf-201803/G3/G3.dREG.peak.full.bed.gz";
file.DREG.G4 <- "../new-rf-201803/G4/G4.dREG.peak.full.bed.gz";
file.DREG.G5 <- "../new-rf-201803/G5/G5.dREG.peak.full.bed.gz";
file.DREG.G6 <- "../new-rf-201803/G6/G6.dREG.peak.full.bed.gz";
file.DREG.G7 <- "../new-rf-201803/G7/G7.dREG.peak.full.bed.gz";
file.DREG.G8 <- "../new-rf-201803/G8/G8.dREG.peak.full.bed.gz";


file.bw.G1 <- c( "../k562/K562_unt.sort.bed.gz_plus.bw", "../k562/K562_unt.sort.bed.gz_minus.bw");
file.bw.G2 <- c( "../k562/groseq_plus.bigWig", "../k562/groseq_minus.bigWig");
file.bw.G3 <- c( "../k562/K562_Nuc_NoRNase_plus.bw", "../k562/K562_Nuc_NoRNase_minus.bw");
file.bw.G4 <- c( "../k562/K562_Nuc_RNase_plus.bw", "../k562/K562_Nuc_RNase_minus.bw");
file.bw.G5 <- c( "../k562/K562_FC_NHS_BRs_normalized_pl.bigWig", "../k562/K562_FC_NHS_BRs_normalized_mn.bigWig");
file.bw.G6 <- c( "../k562/6045_7157_27170_HNHKJBGXX_K562_0min_celastrol10uM_rep1_GB_ATCACG_R1_plus.primary.bw", "../k562/6045_7157_27170_HNHKJBGXX_K562_0min_celastrol10uM_rep1_GB_ATCACG_R1_minus.primary.bw");
file.bw.G7 <- c( "../k562/6045_7157_27176_HNHKJBGXX_K562_0min_celastrol10uM_rep2_GB_CAGATC_R1_plus.primary.bw", "../k562/6045_7157_27176_HNHKJBGXX_K562_0min_celastrol10uM_rep2_GB_CAGATC_R1_minus.primary.bw");
file.bw.G8 <- c( "../k562/SRR182390x_plus.bw", "../k562/SRR182390x_minus.bw");
file.bw.GH <- c( "../HCT116/SRR1105736.7.plus.bw", "../HCT116/SRR1105736.7.minus.bw");
file.bw.GM <- c( "../GM12878/groseq_plus.bigWig", "../GM12878/groseq_minus.bigWig");

options("scipen"=100, "digits"=4)

bed_subtract <- function( score_bed, refer_bed ){
	file.refer <- tempfile(fileext=".bed");
	file.score <- tempfile(fileext=".bed");
	write.table( refer_bed[,c(1:3)], file=file.refer, quote=F, col.names=F, row.names=F, sep="\t");
	write.table( score_bed, file=file.score, quote=F, col.names=F, row.names=F, sep="\t");
	ret <- unique( read.table( pipe(paste("bedtools intersect -a ", file.score,  " -b ", file.refer, " -v ", sep=" ")) ) );
	unlink(file.refer);
	unlink(file.score);
	return(ret);
}


library(bigWig);
pbw <- load.bigWig(file.bw.G1[1]);
mbw <- load.bigWig(file.bw.G1[2]);

## Read in refseq genes.
refGene <- read.table("/fs/cbsudanko/storage/projects/mcf7tamres/annotations/refGene.bed.gz")
refGene <- refGene[grep("random|Un|hap", refGene$V1, invert=TRUE),]
refGene <- refGene[(refGene$V3-refGene$V2)>5000,]

bodies <- refGene
bodies$V2[bodies$V6 == "+"] <-bodies$V2[bodies$V6 == "+"]+1000
bodies$V3[bodies$V6 == "-"] <- bodies$V3[bodies$V6 == "-"]-1000

tb0 <- read.table(file.DREG.G1);
tb0 <- bed_subtract( tb0, bodies );

tb <- tb0[,c(1,6, 6)];
tb[,2] <- tb[,2] -250

LP <- bed.region.bpQuery.bigWig( pbw, tb ) 
LM <- abs(bed.region.bpQuery.bigWig( mbw, tb ));

tb <- tb0[,c(1,6, 6)];
tb[,3] <- tb[,3] + 250
RP <- bed.region.bpQuery.bigWig( pbw, tb ) 
RM <- abs(bed.region.bpQuery.bigWig( mbw, tb ));

if(0)
{
ratio <- (LP+RM+1)/(RP+LM+1)
ratio <- ratio [ratio > 0.025 & ratio < 40]

#hist(log10(ratio), breaks=100);dev.off()
#qqplot( -log10(pnorm( log10(ratio), lower.tail=F, mean=-0.56, sd=0.74)), -log10(runif(NROW(ratio))));segments(0,0,5,5,col="blue");dev.off()
#sd(log10(ratio))
#10^(sd(log10(ratio))*2 - 0.56)
sum(ratio < 1/8.696) 
sum(ratio > 8.696)   
NROW(ratio)          
}

#uni-direction
cat("bi=", NROW( which( RP+LP > 0  & RM+LM > 0 ) ), "fraction=", NROW( which( RP+LP > 0  & RM+LM > 0 ) )/ NROW(LP), "\n" );
# direction
cat("uni=", NROW( which( RP+LP == 0 | RM+LM == 0 ) ), "fraction=", NROW( which( RP+LP ==0  | RM+LM== 0 ) )/ NROW(LP), "\n" );
cat("forward=", NROW( which( RP+LP == 0 ) ), "fraction=", NROW( which( RP+LP ==0  ) )/ NROW(LP), "\n" );
cat("reverse=", NROW( which( RM+LM == 0 ) ), "fraction=", NROW( which( RM+LM ==0  ) )/ NROW(LP), "\n" );


idx <- which( RP+LP > 0  & RM+LM > 0 )
pdf("fig-S5D-uni-direction.pdf")
hist( log((RP+LP)[idx]/(RM+LM)[idx]), breaks=100, main="uni-direction", xlab="log(R/M)");
dev.off();
