
file.G1 <- "../new-model-201712/G1/G1.bedgraph.dREG.peak.score.bed.gz"
file.G2 <- "../new-model-201712/G2/G2.bedgraph.dREG.peak.score.bed.gz"
file.G3 <- "../new-model-201712/G3/G3.bedgraph.dREG.peak.score.bed.gz"
file.G4 <- "../new-model-201712/G4/G4.bedgraph.dREG.peak.score.bed.gz"
file.G5 <- "../new-model-201712/G5/G5.bedgraph.dREG.peak.score.bed.gz"
file.G6 <- "../new-model-201712/G6/G6.bedgraph.dREG.peak.score.bed.gz"
file.G7 <- "../new-model-201712/G7/G7.bedgraph.dREG.peak.score.bed.gz"

system("zcat /fs/cbsudanko/storage/data/hg19/all/ENCODE_tf_peak_calls/wgEncodeRegTfbsClusteredWithCellsV3.bed.gz | grep K562 > K562.tf.bed");

unintersect<-function( file1, file2 )
{
	#Only report those entries in A that have no overlap in B. Restricted by -f and -r.
	tbx <- NULL;
	try( tbx <- read.table(pipe(paste("bedtools intersect -a ", file1, "-b ", file2, " -v"))));
	return(tbx);
}

intersect<-function( file1, file2, flag )
{
	tbx <- NULL;
	try( tbx <- read.table(pipe(paste("bedtools intersect -a ", file1, "-b ", file2, flag ))));
	return(tbx);
}



if(0)
{
	Gr1 <- NROW(unintersect(file.G1, "K562.tf.bed"))/NROW(read.table(file.G1));
	Gr2 <- NROW(unintersect(file.G2, "K562.tf.bed"))/NROW(read.table(file.G2));
	Gr3 <- NROW(unintersect(file.G3, "K562.tf.bed"))/NROW(read.table(file.G3));
	Gr4 <- NROW(unintersect(file.G4, "K562.tf.bed"))/NROW(read.table(file.G4));
	Gr5 <- NROW(unintersect(file.G5, "K562.tf.bed"))/NROW(read.table(file.G5));
	Gr6 <- NROW(unintersect(file.G6, "K562.tf.bed"))/NROW(read.table(file.G6));
	Gr7 <- NROW(unintersect(file.G7, "K562.tf.bed"))/NROW(read.table(file.G7));

	R <- c(Gr1, Gr2, Gr3, Gr4, Gr5, Gr6, Gr7);

	round(R,3)
	#0.166 0.037 0.057 0.087 0.070 0.026 0.027
}

if(1)
{
	novel_TF <- function(file.G)
	{
		file.K562.H3K27ac.peak  <- "/fs/cbsudanko/storage/data/hg19/k562/histones/wgEncodeBroadHistoneK562H3k27acStdAln.bed.gz"
		file.K562.DHS.peak  <- "../k562/K562.bed";

		tb <- unique(read.table(pipe(paste("zcat ", file.K562.H3K27ac.peak , " | bedtools intersect -a" , file.G, "-b - -loj")))[,c(1:5)]);
		idx.novel <- which( unlist(apply(tb, 1, function(x){return(as.character(x[5])==".")})) );
		tb.nov <- tb[ idx.novel, 1:4];

		file1= tempfile();
		write.table(tb.nov, file=file1, quote=F, row.names=F, col.names=F, sep="\t");

		tb <- unique(read.table(pipe(paste("cat ", file.K562.DHS.peak, " | bedtools intersect -a" , file1, "-b - -loj")))[,c(1:5)]);
		idx.novel <- which( unlist(apply(tb, 1, function(x){return(as.character(x[5])==".")})) );
		tb.nov <- tb[ idx.novel, 1:4];

		file1= tempfile();
		write.table(tb.nov, file=file1, quote=F, row.names=F, col.names=F, sep="\t");

		tb.overTF <- unique(intersect( "K562.tf.bed", file1, "-wa" )[,c(1:5)]);
		tb.overPeak <- unique(intersect( file1, "K562.tf.bed", "-wa"))
		
		file1= tempfile();
		write.table(tb.overPeak, file=file1, quote=F, row.names=F, col.names=F, sep="\t");
		file2= tempfile();
		write.table(tb.nov, file=file2, quote=F, row.names=F, col.names=F, sep="\t");
		tb.unoverlap <- unintersect(file2, file1)
		
		show(sort(summary(tb.overTF$V4)));
		#    RDBP    NR2C2     BRF1   POLR3G     BRF2     KAP1   ZNF274    HDAC8     BCL3     ELK1  SMARCB1    HDAC6   GTF2F1
		#       9       15       20       21       32       34       54       57       67       80       81       83       86
		# SMARCA4     BDP1     RFX5     EZH2    GTF2B     CHD1    THAP1   SETDB1   BCLAF1   GTF3C2     NFE2    BACH1    STAT2
		#      90      102      108      136      137      141      156      169      184      184      186      193      193
		#    NRF1    SIRT6     SIX5      SP2      SRF   ZNF263   ZBTB33    SAP30     MXI1    GATA1     TAF1   RPC155    STAT1
		#     208      209      230      230      257      259      263      268      284      295      301      325      325
		#    USF2     TAF7     E2F4    MEF2A SIN3AK20      SP1    KDM5B    RBBP5     ETS1     NFYA     CHD2    HDAC1   ARID3A
		#     332      367      373      378      419      430      491      492      497      498      517      579      612
		#  TRIM28     REST  TBL1XR1   STAT5A      PML    HDAC2    CTCFL     MAFF     MAFK    HMGN3     CBX3    GABPA     NFYB
		#     691      696      729      759      776      777      949      978     1001     1072     1125     1166     1199
		#   NR2F2    FOSL1     UBTF      TBP     SMC3     JUNB     IRF1     PHF8     USF1    CCNT2     ATF1   ZNF143     ATF3
		#    1205     1221     1233     1241     1319     1359     1378     1423     1505     1582     1720     1720     1740
		#  ZBTB7A      JUN      FOS     TAL1      YY1    RAD21  BHLHE40    GATA2     SPI1   POLR2A    CEBPB     E2F6    TEAD4
		#    1744     1775     1815     1850     1855     1881     2039     2083     2098     2337     2408     2600     2667
		#    ELF1    EP300     EGR1     CTCF      MYC      MAZ    RCOR1      MAX     JUND
		#    2698     2882     2952     2974     3224     3256     3576     4766     4941

		return( NROW(tb.overPeak)/NROW(tb.nov) )
	}


	R <- c( novel_TF(file.G1), novel_TF(file.G2), novel_TF(file.G3), novel_TF(file.G4), novel_TF(file.G5), novel_TF(file.G6), novel_TF(file.G7) )
}



