options("scipen"=100, "digits"=4)
options(width=250);

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

bed_intersect <- function( score_bed, refer_bed ){
	file.refer <- tempfile(fileext=".bed");
	file.score <- tempfile(fileext=".bed");
	write.table( refer_bed[,c(1:3)], file=file.refer, quote=F, col.names=F, row.names=F, sep="\t");
	write.table( score_bed, file=file.score, quote=F, col.names=F, row.names=F, sep="\t");
	ret <- unique( read.table( pipe(paste("bedtools intersect -a ", file.score,  " -b ", file.refer, " -wa ", sep=" ")) ) );
	unlink(file.refer);
	unlink(file.score);
	return(ret);
}

bed_coverage <- function( score_bed, refer_bed ){
	file.refer <- tempfile(fileext=".bed");
	file.score <- tempfile(fileext=".bed");
	write.table( refer_bed[,c(1:3)], file=file.refer, quote=F, col.names=F, row.names=F, sep="\t");
	write.table( score_bed, file=file.score, quote=F, col.names=F, row.names=F, sep="\t");

	ret <- read.table( pipe(paste("bedtools coverage -counts -a ", file.score,  " -b ", file.refer, sep=" ")));

	unlink(file.refer);
	unlink(file.score);
	return(ret);
}


file.black.bed <- "/fs/cbsudanko/storage/data/hg19/all/encode_blacklist/hg19.blacklist.bed";

file.DREG.G1 <- "../new-rf-201803/G1/G1.dREG.peak.full.bed.gz";
file.DREG.G2 <- "../new-rf-201803/G2/G2.dREG.peak.full.bed.gz";
file.DREG.G3 <- "../new-rf-201803/G3/G3.dREG.peak.full.bed.gz";
file.DREG.G4 <- "../new-rf-201803/G4/G4.dREG.peak.full.bed.gz";
file.DREG.G5 <- "../new-rf-201803/G5/G5.dREG.peak.full.bed.gz";
file.DREG.G6 <- "../new-rf-201803/G6/G6.dREG.peak.full.bed.gz";
file.DREG.G7 <- "../new-rf-201803/G7/G7.dREG.peak.full.bed.gz";
file.DREG.G8 <- "../new-rf-201803/G8/G8.dREG.peak.full.bed.gz";


file.dREG <- c(file.DREG.G1, file.DREG.G2, file.DREG.G3, file.DREG.G4, file.DREG.G5, file.DREG.G6, file.DREG.G7, file.DREG.G8)

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


## BroadHmm has too noises.
if(0)
{
	# V4
	#  1_Active_Promoter 
	#  2_Weak_Promoter 
	#  3_Poised_Promoter
	#  4_Strong_Enhancer 
	#  5_Strong_Enhancer
	#  6_Weak_Enhancer   
	#  7_Weak_Enhancer   
	#  8_Insulator       
	#  9_Txn_Transition  
	#  10_Txn_Elongation
	#  11_Weak_Txn      
	#  12_Repressed      
	#  13_Heterochrom/lo 
	#  14_Repetitive/CNV 
	#  15_Repetitive/CNV 

	#1st try
	#tb.hetero.K562 <- read.table(pipe("zcat wgEncodeBroadHmmK562HMM.bed.gz       | grep '13_' | bedtools merge -i -"));
	#tb.hetero.GM <- read.table(pipe("zcat wgEncodeBroadHmmGm12878HMM.bed.gz | grep '13_' | bedtools merge -i -"));

	##NO ChromHMM file for HCT116
	##tb.hetero.HCT <- read.table(pipe("zcat wgEncodeBroadHmmK562HMM.bed.gz   | grep '13_' | bedtools merge -i -"));

	#2nd try
	#tb.hetero.K562 <- read.table(pipe("zcat wgEncodeBroadHmmK562HMM.bed.gz       | grep '12_' | bedtools merge -i -"));
	#tb.hetero.GM <- read.table(pipe("zcat wgEncodeBroadHmmGm12878HMM.bed.gz | grep '12_' | bedtools merge -i -"));

	##NO ChromHMM file for HCT116
	##tb.hetero.HCT <- read.table(pipe("zcat wgEncodeBroadHmmK562HMM.bed.gz   | grep '12_' | bedtools merge -i -"));

	#3rd try
	tb.hetero.K562 <- read.table(pipe("zcat wgEncodeBroadHmmK562HMM.bed.gz  | grep -E '12_|13_' | bedtools merge -i -"));
	tb.hetero.GM <- read.table(pipe("zcat wgEncodeBroadHmmGm12878HMM.bed.gz | grep -E '12_|13_' | bedtools merge -i -"));

	##NO ChromHMM file for HCT116
	#tb.hetero.HCT <- read.table(pipe("zcat wgEncodeBroadHmmK562HMM.bed.gz   | grep -E '12_|13_' | bedtools merge -i -"));
}

tb.hetero.K562 <- read.table(pipe("zcat /fs/cbsudanko/storage/data/hg19/k562/histones/wgEncodeBroadHistoneK562H3k27me3StdPk.broadPeak.gz"));
tb.hetero.GM <- read.table(pipe("zcat /fs/cbsudanko/storage/data/hg19/gm12878/histones/wgEncodeBroadHistoneGm12878H3k27me3StdPkV2.broadPeak.gz"));

dREG.intersect.hetero <- list();

df.hetero<-c();
tb.refer <- read.table(file.K562.DHS.peak);
for(fdREG in file.dREG)
{
	tb.dREG <- read.table( pipe(paste("bedtools subtract -A -a ", fdREG, "-b", file.black.bed)));

	tb.DHSminus <- bed_subtract(tb.dREG, tb.refer);
	tb.overlap <- bed_intersect( tb.DHSminus, tb.hetero.K562);
	tb.overlap2 <- bed_intersect( tb.hetero.K562, tb.DHSminus );

	tba.overlap <- bed_intersect( tb.dREG, tb.hetero.K562 );
	tba.overlap2 <- bed_intersect( tb.hetero.K562, tb.dREG );
    

    df.hetero<-rbind(df.hetero, c( NROW(tb.overlap), NROW(tb.overlap)/ NROW(tb.dREG),  NROW(tb.overlap2),  NROW(tb.overlap2)/NROW(tb.hetero.K562), 
                                   NROW(tba.overlap), NROW(tba.overlap)/ NROW(tb.dREG),  NROW(tba.overlap2), NROW(tba.overlap2)/NROW(tb.hetero.K562) ) );

    dREG.intersect.hetero[[fdREG]] <- tba.overlap
	print(fdREG);
	cat ( NROW(tb.overlap), NROW(tb.overlap)/ NROW(tb.dREG),  NROW(tb.overlap2)/NROW(tb.hetero.K562), "\n" );
}


unlink(file.K562.DHS.peak);

#-----------------
# GM12878
#-----------------

file.GM.H3K27ac.peak  <- "/fs/cbsudanko/storage/data/hg19/gm12878/histones/wgEncodeBroadHistoneGm12878H3k27acStdPk.broadPeak.gz"
file.GM12878.DHS.peak <- "./uw.merge.narrowPeak.bed"

make_GM12878_bed<-function( file.nonneg,dist=100 )
{
	nonneg_bed <- read.table("/fs/cbsudanko/storage/data/hg19/gm12878/dnase/uw.merge.narrowPeak.bed");
	nonneg_bed[,2] <- nonneg_bed[,2] - dist
	idx.mis <- which(nonneg_bed[,2]<0);
	if(length(idx.mis)>0) nonneg_bed[idx.mis,2] <- 0;
	nonneg_bed[,3] <- nonneg_bed[,3] + dist
	write.table( nonneg_bed, file=file.nonneg, quote=F, row.name=F, col.names=F, sep="\t" );
}

make_GM12878_bed( file.GM12878.DHS.peak );
tb.refer <- read.table(file.GM12878.DHS.peak);


if(1)
{
	file.DREG.GM <- "../new-rf-201803/GM/GM.dREG.peak.full.bed.gz";
	tb.dREG <- read.table( pipe(paste("bedtools subtract -A -a ", file.DREG.GM, "-b", file.black.bed)));

	tb.DHSminus <- bed_subtract(tb.dREG, tb.refer);
	tb.overlap <- bed_intersect( tb.DHSminus, tb.hetero.GM );
	tb.overlap2 <- bed_intersect( tb.hetero.GM, tb.DHSminus );

	tba.overlap <- bed_intersect(  tb.dREG, tb.hetero.GM );
	tba.overlap2 <- bed_intersect( tb.hetero.GM, tb.dREG );
	df.hetero<-rbind(df.hetero, c( NROW(tb.overlap), NROW(tb.overlap)/ NROW(tb.dREG),  NROW(tb.overlap2),  NROW(tb.overlap2)/NROW(tb.hetero.GM), 
									   NROW(tba.overlap), NROW(tba.overlap)/ NROW(tb.dREG),  NROW(tba.overlap2), NROW(tba.overlap2)/NROW(tb.hetero.GM) ) );
	dREG.intersect.hetero[[file.DREG.GM]] <- tba.overlap;

	print(file.DREG.GM);
	cat( NROW(tb.overlap), NROW(tb.overlap)/ NROW(tb.dREG),  NROW(tb.overlap2)/NROW(tb.hetero.GM), "\n" );
}

unlink(file.GM12878.DHS.peak);

file.dREG <-  c(file.dREG, file.DREG.GM);

#-----------------
# HCT116
#-----------------

if(0)
{
file.HCT116.DHS.peak <- "./HCT116.narrowPeak.bed"

make_HCT116_bed<-function( file.nonneg,dist=100 )
{
	nonneg_bed <- read.table("../HCT116/hg19.dNase.HCT116.broad.peak.bed");
	nonneg_bed[,2] <- nonneg_bed[,2] - dist
	idx.mis <- which(nonneg_bed[,2]<0);
	if(length(idx.mis)>0) nonneg_bed[idx.mis,2] <- 0;
	nonneg_bed[,3] <- nonneg_bed[,3] + dist
	write.table( nonneg_bed, file=file.nonneg, quote=F, row.name=F, col.names=F, sep="\t" );
}

make_HCT116_bed( file.HCT116.DHS.peak );
tb.refer <- read.table(file.HCT116.DHS.peak);
file.DREG.GH <- "../new-rf-201803/GH/GH.dREG.peak.full.bed.gz";

	tb.dREG <- read.table( pipe(paste("bedtools subtract -A -a ", file.DREG.GH, "-b", file.black.bed)));
	tb.DHSminus <- bed_subtract( tb.dREG, tb.refer );
	tb.overlap <- bed_intersect( tb.DHSminus, tb.hetero.HCT );
	tb.overlap2 <- bed_intersect( tb.hetero.HCT , tb.DHSminus );

	tba.overlap <- bed_intersect( tb.dREG, tb.hetero.HCT );
	tba.overlap2 <- bed_intersect( tb.hetero.HCT, tb.dREG );
    df.hetero<-rbind(df.hetero, c( NROW(tb.overlap), NROW(tb.overlap)/ NROW(tb.dREG),  NROW(tb.overlap2),  NROW(tb.overlap2)/NROW(tb.hetero.HCT), 
                                   NROW(tba.overlap), NROW(tba.overlap)/ NROW(tb.dREG),  NROW(tba.overlap2), NROW(tba.overlap2)/NROW(tb.hetero.HCT) ) );
	
	print(file.DREG.GH);
	cat( NROW(tb.overlap), NROW(tb.overlap)/ NROW(tb.dREG),  NROW(tb.overlap2)/NROW(tb.hetero.HCT), "\n" );


unlink(file.HCT116.DHS.peak);
}


show(df.hetero);

#use BraodHmm
if(0)
{
# 13_Heterochrom/lo total:78279
#      [,1]     [,2] [,3]     [,4]  [,5]    [,6] [,7]     [,8]
# [G1,] 2880 0.038149 2500 0.031937  4923 0.06521 4310 0.055059
# [G2,]  143 0.005484  140 0.001788   300 0.01150  291 0.003717
# [G3,]  675 0.012288  641 0.008189  1523 0.02773 1435 0.018332
# [G4,]  523 0.010437  465 0.005940  1208 0.02411 1092 0.013950
# [G5,]  602 0.012857  550 0.007026  1368 0.02922 1250 0.015969
# [G6,]  156 0.004446  153 0.001955   428 0.01220  415 0.005302
# [G7,]  146 0.004216  143 0.001827   409 0.01181  397 0.005072
# [G8,] 1402 0.032537 1295 0.016543  2092 0.04855 1951 0.024924
# [GM,] 1816 0.025543 1624 0.021621  2892 0.04068 2607 0.034708


# 12_Repressed total:38902
#      [,1]     [,2] [,3]      [,4] [,5]     [,6] [,7]     [,8]
# [G1,]  688 0.009113  634 0.0162974 1514 0.020055 1379 0.035448
# [G2,]   46 0.001764   45 0.0011568  101 0.003873  100 0.002571
# [G3,]  126 0.002294  123 0.0031618  508 0.009248  493 0.012673
# [G4,]  122 0.002435  115 0.0029561  445 0.008880  412 0.010591
# [G5,]  130 0.002776  126 0.0032389  466 0.009952  448 0.011516
# [G6,]   39 0.001111   37 0.0009511  163 0.004645  160 0.004113
# [G7,]   44 0.001271   44 0.0011310  159 0.004591  158 0.004061
# [G8,]  451 0.010466  413 0.0106164  699 0.016222  648 0.016657
# [GM,]  270 0.003798  262 0.0102814  431 0.006062  416 0.016325


# 12_Repressed + 13_Heterochrom/lo total: 73338
#      [,1]     [,2] [,3]     [,4]  [,5]    [,6]  [,7]     [,8]
# [G1,] 3532 0.046785 3062 0.041752  6375 0.08444  5546 0.075622
# [G2,]  189 0.007247  184 0.002509   397 0.01522   386 0.005263
# [G3,]  792 0.014418  755 0.010295  2008 0.03655  1896 0.025853
# [G4,]  639 0.012751  572 0.007800  1638 0.03269  1485 0.020249
# [G5,]  725 0.015484  669 0.009122  1817 0.03880  1677 0.022867
# [G6,]  194 0.005528  188 0.002563   587 0.01673   570 0.007772
# [G7,]  188 0.005429  184 0.002509   564 0.01629   551 0.007513
# [G8,] 1824 0.042330 1674 0.022826  2746 0.06373  2541 0.034648
# [GM,] 2066 0.029059 1859 0.027667  3301 0.04643  2975 0.044276
}

#      [,1]     [,2] [,3]     [,4]  [,5]   [,6] [,7]    [,8]
# [1,] 2315 0.030665 1610 0.018281 16293 0.2158 5321 0.06042
# [2,]  210 0.008053  195 0.002214  4824 0.1850 2309 0.02622
# [3,]  535 0.009739  452 0.005132 11172 0.2034 3874 0.04399
# [4,]  589 0.011754  475 0.005393  9563 0.1908 3352 0.03806
# [5,]  534 0.011404  474 0.005382  9614 0.2053 3671 0.04168
# [6,]  185 0.005272  180 0.002044  6578 0.1875 2682 0.03045
# [7,]  187 0.005400  191 0.002169  6469 0.1868 2621 0.02976
# [8,] 1215 0.028197  944 0.010719  8119 0.1884 3316 0.03765
# [9,] 4686 0.065910 2168 0.052286 18933 0.2663 4339 0.10464


POLR2A.bed <- read.table(pipe("zcat /fs/cbsudanko/storage/data/hg19/all/ENCODE_tf_peak_calls/wgEncodeRegTfbsClusteredWithCellsV3.bed.gz | grep POLR2A - | grep  K562 "))
POLR2A.bed <- POLR2A.bed[,c(1:3)]
RPC155.bed <- read.table(pipe("zcat /fs/cbsudanko/storage/data/hg19/all/ENCODE_tf_peak_calls/wgEncodeRegTfbsClusteredWithCellsV3.bed.gz | grep RPC155 - | grep  K562 "))
RPC155.bed <- RPC155.bed[,c(1:3)]

df.dREG.hetero=array(0, dim=c(NROW(dREG.intersect.hetero),8));

for(i in 1:NROW(dREG.intersect.hetero))
{
	df.dREG.hetero[i,1] <- NROW( read.table( pipe(paste("bedtools subtract -A -a ", file.dREG[i], "-b", file.black.bed))))
	df.dREG.hetero[i,2] <- NROW( dREG.intersect.hetero[[i]] );
	df.dREG.hetero[i,3] <- NROW( POLR2A.bed );
	df.dREG.hetero[i,4] <- NROW( RPC155.bed );
	df.dREG.hetero[i,5] <- NROW( bed_intersect( read.table( pipe(paste("bedtools subtract -A -a ", file.dREG[i], "-b", file.black.bed))) ,POLR2A.bed));
	df.dREG.hetero[i,6] <- NROW( bed_intersect( read.table( pipe(paste("bedtools subtract -A -a ", file.dREG[i], "-b", file.black.bed))), RPC155.bed));
	df.dREG.hetero[i,7] <- NROW( bed_intersect( dREG.intersect.hetero[[i]], POLR2A.bed));
	df.dREG.hetero[i,8] <- NROW( bed_intersect( dREG.intersect.hetero[[i]], RPC155.bed));
}


colnames(df.dREG.hetero) <- c("Total.dREG", "Total.Hetero", "Total.POLR2A", "Total.RPC155", "dREG*POLR2A", "dREG*RPC155", "Hetero*POLR2A", "Hetero*RPC155" )
rownames(df.dREG.hetero) <- c("G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8", "GM")
show(df.dREG.hetero);

if(1)
{
library(scales)
	tb <- bed_coverage(tb.hetero.K562, read.table( pipe(paste("bedtools subtract -A -a ", file.dREG[1], "-b", file.black.bed)))  );
	tb2 <- data.frame(tb$V3-tb$V2, tb$V10)
	colnames(tb2) <-c("width", "TIR")
	loessMod10 <- loess(TIR~width, data=tb2, span=0.1) # 10% smoothing span
    tb2<- tb2[tb2$width>=10000,]
    
#	pdf("fig-4B-logx.pdf")     
#	plot(tb2[,1], tb2[,2], xlab="Length of Heterochromatin Peak", ylab="TIR count", col=alpha("#99a1af", 0.4), log="x", pch=20, ylim=c(0,20), xlim=c(10000, 500000) )
#	smoothed10 <- predict(loessMod10, seq(10000, max(tb2[,1]), 1000) )
#	lines(seq(10000, max(tb2[,1]), 1000), smoothed10 );
#	dev.off()

library(ggplot2); 
library(data.table); 
library(parallel)

    x <- seq( floor(log(10000)), ceiling(log(500000)), 0.01 )
	tb.ci <- rbindlist(mclapply(x , function(i){
	    TIR <- tb2[tb2$width>=exp(i) & tb2$width < exp(i+0.01), "TIR" ]
		return(  data.frame(x=exp(i), c05=quantile(TIR, 0.05), mean=mean(TIR), c95=quantile(TIR, 0.95) ) );
	}, mc.cores=1));
	tb.ci <- tb.ci[!is.na(tb.ci$c05) & !is.na(tb.ci$c95),]
    
    tb20 <- smooth.spline(x=tb.ci$x, y=tb.ci$c05)
    tb21 <- smooth.spline(x=tb.ci$x, y=tb.ci$c95)
    tb22 <- smooth.spline(x=tb.ci$x, y=tb.ci$mean)
	tb2.ci <- data.frame(x=tb20$x, c25=tb20$y, c75=tb21$y, y=tb22$y);

    ggplot(data=tb2.ci, aes(x, y), log="x" ) +
        ylab("TIR count") +
        xlab("Length of Heterochromatin Peak") + 
        scale_x_continuous(trans="log", limits=c(10000, 500000), breaks=c(10000, 20000, 50000, 100000, 250000, 500000) )+
        ylim(-1, 20)+
        geom_ribbon( aes(ymin=c25, ymax=c75), linetype=1, fill = "#5070fb", alpha=0.5, size=0.1)+
        geom_line() + 
        theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))
        
    ggsave("fig-4B-logx.pdf", device="pdf", width=4, height=4)

}

if(1)
{
	tb.TIR <- bed_intersect(read.table( pipe(paste("bedtools subtract -A -a ", file.dREG[1], "-b", file.black.bed))), tb.hetero.K562 );
	tb.non <- bed_subtract(read.table( pipe(paste("bedtools subtract -A -a ", file.dREG[1], "-b", file.black.bed))), tb.TIR );

	file.bw.G1 <- c( "../k562/K562_unt.sort.bed.gz_plus.bw", "../k562/K562_unt.sort.bed.gz_minus.bw");
	file.G1.plus <- file.bw.G1[1];
	file.G1.minus <- file.bw.G1[2];

	tb.TIR[,3] <- tb.TIR[,3] + 250
	tb.TIR[,2] <- tb.TIR[,2] - 250
	tb.non[,3] <- tb.non[,3] + 250
	tb.non[,2] <- tb.non[,2] - 250
	
	library(bigWig);
	bw.plus <- load.bigWig( file.G1.plus);
	bw.minus <- load.bigWig( file.G1.minus);
	reads.TIR <- bed.region.bpQuery.bigWig( bw.plus, tb.TIR[,c(1:3)], op = "sum", abs.value = TRUE) + 
	             bed.region.bpQuery.bigWig( bw.minus,tb.TIR[,c(1:3)], op = "sum", abs.value = TRUE);
	reads.non <- bed.region.bpQuery.bigWig( bw.plus, tb.non[,c(1:3)], op = "sum", abs.value = TRUE) + 
	             bed.region.bpQuery.bigWig( bw.minus,tb.non[,c(1:3)], op = "sum", abs.value = TRUE);
	unload.bigWig(bw.plus);
	unload.bigWig(bw.minus);
	reads <- data.frame(count=c(reads.TIR, reads.non)+1, type=c(rep("TIR in Heterochromatin", NROW(reads.TIR)), rep("TIR not in Heterochromatin", NROW(reads.non)) ) );

	pdf("fig-4C.pdf")
	boxplot(count~type,data=reads, main="",   xlab="", ylab="Read counts", log="y" ) 
	dev.off();
	
}


#test POLR2A is more frequent in heterochromatin in G1
x <- array(df.dREG.hetero[1,c("Total.dREG","Total.Hetero", "dREG*POLR2A", "Hetero*POLR2A")], dim=c(2,2))
x[1,] <- x[1,] - x[2,]
fisher.test(x)

#test POLR2A is more frequent in heterochromatin in GM
x <- array(df.dREG.hetero[9,c("Total.dREG","Total.Hetero", "dREG*POLR2A", "Hetero*POLR2A")], dim=c(2,2))
x[1,] <- x[1,] - x[2,]
fisher.test(x)

#test RPC155 is more frequent in heterochromatin in G1
x <- array(df.dREG.hetero[1,c("Total.dREG","Total.Hetero", "dREG*RPC155", "Hetero*RPC155")], dim=c(2,2))
x[1,] <- x[1,] - x[2,]
fisher.test(x)

#test RPC155 is more frequent in heterochromatin in GM
x <- array(df.dREG.hetero[9,c("Total.dREG","Total.Hetero", "dREG*RPC155", "Hetero*RPC155")], dim=c(2,2))
x[1,] <- x[1,] - x[2,]
fisher.test(x)
