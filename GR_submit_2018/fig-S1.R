## get correlations between MCF-7 samples

## Read in refseq genes.
refGene <- read.table("/fs/cbsudanko/storage/projects/mcf7tamres/annotations/refGene.bed.gz")
refGene <- refGene[grep("random|Un|hap", refGene$V1, invert=TRUE),]
refGene <- refGene[(refGene$V3-refGene$V2)>5000,]

bodies <- refGene
bodies$V2[bodies$V6 == "+"] <-bodies$V2[bodies$V6 == "+"]+1000
bodies$V3[bodies$V6 == "-"] <- bodies$V3[bodies$V6 == "-"]-1000

tss <- refGene
tss$V3[tss$V6 == "+"] <-tss$V2[tss$V6 == "+"]+1000
tss$V2[tss$V6 == "-"] <- tss$V3[tss$V6 == "-"]-1000

## Read in bigWigs.
require(bigWig)

G1_pl <- load.bigWig("../k562/K562_unt.sort.bed.gz_plus.bw")
G1_mn <- load.bigWig("../k562/K562_unt.sort.bed.gz_minus.bw")

G2_pl <- load.bigWig("../k562/groseq_plus.bigWig")
G2_mn <- load.bigWig("../k562/groseq_minus.bigWig")

G3_pl <- load.bigWig("../k562/K562_Nuc_NoRNase_plus.bw")
G3_mn <- load.bigWig("../k562/K562_Nuc_NoRNase_minus.bw")

G4_pl <- load.bigWig("../k562/K562_Nuc_RNase_plus.bw")
G4_mn <- load.bigWig("../k562/K562_Nuc_RNase_minus.bw")

G5_pl <- load.bigWig("../k562/K562_FC_NHS_BRs_normalized_pl.bigWig")
G5_mn <- load.bigWig("../k562/K562_FC_NHS_BRs_normalized_mn.bigWig")

G6_pl <- load.bigWig("../k562/6045_7157_27170_HNHKJBGXX_K562_0min_celastrol10uM_rep1_GB_ATCACG_R1_plus.primary.bw");
G6_mn <- load.bigWig("../k562/6045_7157_27170_HNHKJBGXX_K562_0min_celastrol10uM_rep1_GB_ATCACG_R1_minus.primary.bw");

G7_pl <- load.bigWig("../k562/6045_7157_27176_HNHKJBGXX_K562_0min_celastrol10uM_rep2_GB_CAGATC_R1_plus.primary.bw");
G7_mn <- load.bigWig("../k562/6045_7157_27176_HNHKJBGXX_K562_0min_celastrol10uM_rep2_GB_CAGATC_R1_minus.primary.bw");

## Count reads in each ...
G1 <- bed6.region.bpQuery.bigWig(G1_pl, G1_mn, bodies, abs.value = TRUE)/ (bodies$V3-bodies$V2) * 1000/ (G1_pl$basesCovered * G1_pl$mean + abs(G1_mn$basesCovered * G1_mn$mean)) * 1e6
G2 <- bed6.region.bpQuery.bigWig(G2_pl, G2_mn, bodies, abs.value = TRUE)/ (bodies$V3-bodies$V2) * 1000/ (G2_pl$basesCovered * G2_pl$mean + abs(G2_mn$basesCovered * G2_mn$mean)) * 1e6
G3 <- bed6.region.bpQuery.bigWig(G3_pl, G3_mn, bodies, abs.value = TRUE)/ (bodies$V3-bodies$V2) * 1000/ (G3_pl$basesCovered * G3_pl$mean + abs(G3_mn$basesCovered * G3_mn$mean)) * 1e6
G4 <- bed6.region.bpQuery.bigWig(G4_pl, G4_mn, bodies, abs.value = TRUE)/ (bodies$V3-bodies$V2) * 1000/ (G4_pl$basesCovered * G4_pl$mean + abs(G4_mn$basesCovered * G4_mn$mean)) * 1e6
G5 <- bed6.region.bpQuery.bigWig(G5_pl, G5_mn, bodies, abs.value = TRUE)/ (bodies$V3-bodies$V2) * 1000/ (G5_pl$basesCovered * G5_pl$mean + abs(G5_mn$basesCovered * G5_mn$mean)) * 1e6
G6 <- bed6.region.bpQuery.bigWig(G6_pl, G6_mn, bodies, abs.value = TRUE)/ (bodies$V3-bodies$V2) * 1000/ (G6_pl$basesCovered * G6_pl$mean + abs(G6_mn$basesCovered * G6_mn$mean)) * 1e6
G7 <- bed6.region.bpQuery.bigWig(G7_pl, G7_mn, bodies, abs.value = TRUE)/ (bodies$V3-bodies$V2) * 1000/ (G7_pl$basesCovered * G7_pl$mean + abs(G7_mn$basesCovered * G7_mn$mean)) * 1e6

gene_body_counts <- cbind(G1, G2, G3, G4, G5, G6, G7 )

## Make some plots
plot(G1, G2)
pairs(gene_body_counts, log="xy")
cor(gene_body_counts, method="spearman")
cor(gene_body_counts, method="spearman")

write.csv(cbind(refGene, gene_body_counts), "gene_body_counts.csv")

## Cluster changed genes...
yb.sig.pal <- function(n, scale=10) {
 ints<- c(0:(n-1))/(n-1)   ## Linear scale from 0:1 x N values.
 ints<- 1/(1+exp(scale*(0.5-ints)))## Transfer to sigmoidal scale.
 b<- min(ints)
 m<- 2*b/(n-1)
 ints<- ints+(m*c(0:(n-1)) -b)## Transform by linear function to fill colors out to maxes.

 ## Transfer to colorspace.
 # Yellow: 255, 255, 0
 # White:  255, 255, 255
 # Blue:   0, 0, 255
 YW <- ints[ints < 0.5] *2
 WB <- (ints[ints >= 0.5]-0.5) *2
 YW[YW<0] <- 0; WB[WB>1] <- 1
 c(rgb(1, 1, YW), rgb(1-WB, 1-WB, 1))
}

drawCor <- function(indx) {
library(lattice);
library(latticeExtra);

	rpkm_df <- as.matrix(gene_body_counts[,indx]) # as.matrix(ca[,indx])#/(ca[,"mapSize"]) ## "Good?!"  Remove H2-U, H3-PI, C2-U+PI, M1-PI

	cond <- c(1,1,1,1,2,2,3,4,5,6,7,3,4,5,6,7)[indx]#"", "", "", "", "", "")# Condition[indx]
	spec <- c(1,1,2,2,1,2,1,1,1,1,1,2,2,2,2,2)[indx]#"", "", "", "", "", "")
	labs <- colnames(gene_body_counts)[indx]

	cc <- cor(rpkm_df, method="spearman")
	#clu <- agnes(t(rpkm_df))

	pal3 <- c("#E03CE9", "#17B92B", "#E6350D", "#6FD2F0", "#F9F77F", "#5B6C0C", "#68003D", "#310F08")

	## Print dendrogram and heatmap with latticeExtra.
	# hc1 <- agnes(1-cc, diss=TRUE, method="ward")
	# hc1 <- hclust(dist(t(rpkm_df), method = "canberra"))
	 hc1 <- hclust(dist(cc, method = "euclidean"),  method="single")
	 hc1 <- as.dendrogram(hc1)
	 ord.hc1 <- order.dendrogram(hc1)
	 hc2 <- reorder(hc1, cond[ord.hc1])
	 ord.hc2 <- order.dendrogram(hc2)
	 #region.colors <- trellis.par.get("superpose.polygon")$col

	 pl <- levelplot((cc)[ord.hc2, ord.hc2], col.regions= yb.sig.pal(100, scale=3), xlab="", ylab="", #rev(cm.colors(100)),  # #c("white", "yellow", "blue") # c("#E9F231", "#B1EC2C", "#5DBDEF")
		 colorkey = list(space="left", labels=list(cex=1.5)),
		 scales = list(x= list(rot=90, cex=1.5, labels=labs[ord.hc2]), y=list(draw=FALSE)), #scales = list(x = list(rot = 90)),
		 legend = list(
			right = list(fun = dendrogramGrob,
				 args = list(x = hc2, ord = ord.hc2, side = "right", #lwd=2,
				 size = 7, size.add = 0.5,
				 add = list(rect = list(col = "transparent", fill = pal3[c(1:8)][cond])),
				 type = "rectangle")),
			top = list(fun = dendrogramGrob,
				 args = list(x = hc2, ord = ord.hc2, side = "top", #lwd=2,
				 size = 1, size.add = 0.5,
				 add = list(rect = list(col = "transparent", fill = pal3[2:8][spec])),
				 type = "rectangle"))
				 ))
	 print(pl)
}

if(0)
{
	pdf("fig-S1-cor.pdf")
		drawCor(c(1:3,5,6))
	dev.off()
}

drawPairs <- function(indx,log=FALSE) {
	rpkm <- as.matrix(gene_body_counts[,indx]) # as.matrix(ca[,indx])#/(ca[,"mapSize"]) ## "Good?!"  Remove H2-U, H3-PI, C2-U+PI, M1-PI

	labs <- colnames(gene_body_counts)[indx]

	if (log)
		pairs(log(rpkm[,c(5,2,4,1,3)]), pch=21, horInd=c(1:5), verInd=c(5,4,3,2,1))
	else
		pairs(rpkm[,c(5,2,4,1,3)], pch=21, horInd=c(1:5), verInd=c(5,4,3,2,1))

    library(lattice)
	if (log)
		splom(~(log(rpkm[,c(5,2,4,1,3)])))
	else
		splom(~rpkm[,c(5,2,4,1,3)]);
}


if(0)
{
	pdf("fig-S1-pairs.pdf")
		drawPairs(c(1:3,5,6), log=TRUE)
	dev.off()
}

drawCorWithPairs <- function(indx, log=FALSE) {
	library(lattice);
	library(grid);
	library(latticeExtra);
	source("levelplot2.R");

assignInNamespace("panel.levelplot", panel.levelplot, ns="lattice")


	rpkm_df <- as.matrix(gene_body_counts[,indx]) # as.matrix(ca[,indx])#/(ca[,"mapSize"]) ## "Good?!"  Remove H2-U, H3-PI, C2-U+PI, M1-PI

	cond <- c(1,1,1,1,2,2,3,4,5,6,7,3,4,5,6,7)[indx]#"", "", "", "", "", "")# Condition[indx]
	spec <- c(1,1,2,2,1,2,1,1,1,1,1,2,2,2,2,2)[indx]#"", "", "", "", "", "")
	labs <- colnames(gene_body_counts)[indx]

	cc <- cor(rpkm_df, method="spearman")
	#clu <- agnes(t(rpkm_df))

	#pal3 <- c("#E03CE9", "#17B92B", "#E6350D", "#6FD2F0", "#F9F77F", "#5B6C0C", "#68003D", "#310F08")
	pal3 <- c("#A0A0A0", "#A0A0A0", "#A0A0A0", "#A0A0A0", "#A0A0A0", "#A0A0A0", "#A0A0A0", "#A0A0A0")

	## Print dendrogram and heatmap with latticeExtra.
	# hc1 <- agnes(1-cc, diss=TRUE, method="ward")
	# hc1 <- hclust(dist(t(rpkm_df), method = "canberra"))
	 hc1 <- hclust(dist(cc, method = "euclidean"),  method="single")
	 hc1 <- as.dendrogram(hc1)
	 ord.hc1 <- order.dendrogram(hc1)
	 hc2 <- reorder(hc1, cond[ord.hc1])
	 ord.hc2 <- order.dendrogram(hc2)
	 #region.colors <- trellis.par.get("superpose.polygon")$col

	 g_pairs <<- if (log) log(rpkm_df[,ord.hc2]+0) else rpkm_df[,ord.hc2];
	 #g_pairs_col <<- "#0080FF";
	 g_pairs_col <<- "#808080";
	 g_pairs_size <<- 0.02;

	 pl <- lattice:::levelplot((cc)[ord.hc2, ord.hc2], xlab="", ylab="",
	     col.regions= yb.sig.pal(100, scale=3),
		 colorkey = list(space="left", labels=list(cex=1.5)),
		 scales = list(x= list(rot=90, cex=1.5, labels=labs[ord.hc2]), y=list(draw=FALSE)),
		 legend = list(
			right = list(fun = dendrogramGrob,
				 args = list(x = hc2, ord = ord.hc2, side = "right", #lwd=2,
				 size = 7, size.add = 0.5,
				 add = list(rect = list(col = "black", fill = pal3[c(1:8)][cond])), #col = "transparent"
				 type = "rectangle")),
			top = list(fun = dendrogramGrob,
				 args = list(x = hc2, ord = ord.hc2, side = "top", #lwd=2,
				 size = 1, size.add = 0.5,
				 add = list(rect = list(col = "black", fill = pal3[2:8][spec])),
				 type = "rectangle"))
				 ))
	 print(pl)
}

if(1)
{
	source("levelplot2.R")
	pdf("fig-S1-corpairs.pdf")
		drawCorWithPairs(c(1:3,5,6), log=TRUE)
	dev.off()
}

unlink("gene_body_counts.csv");
