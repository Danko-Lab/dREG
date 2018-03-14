file.H3k27ac  <- "/fs/cbsudanko/storage/data/hg19/k562/histones/wgEncodeBroadHistoneK562H3k27acStdAln.bed.gz"
dhs_folder <- "/workdir/tc532/manuscript_code/database/merged/min150_sorted_data"

unintersect<-function( bedorfile1, bedorfile2 )
{
	bedfile1 <- bedorfile1;
	bedfile2 <- bedorfile2;

	if(class(bedorfile1)=="data.frame")
	{
		bedfile1 = tempfile();
		write.table(bedorfile1, file=bedfile1, quote=F, col.names=F, row.names=F, sep="\t");
	}

	if(class(bedorfile2)=="data.frame")
	{
		bedfile2 = tempfile();
		write.table(bedorfile2, file=bedfile2, quote=F, col.names=F, row.names=F, sep="\t");
	}

	#Only report those entries in A that have no overlap in B. Restricted by -f and -r.
	tbx <- NULL;
	try( tbx <- read.table(pipe(paste("bedtools intersect -a ", bedfile1, "-b ", bedfile2, " -v"))));

	if(class(bedorfile1)=="data.frame") unlink(bedfile1);
	if(class(bedorfile2)=="data.frame") unlink(bedfile2);

	return(tbx);
}

intersect<-function( bedorfile1, bedorfile2 )
{
	bedfile1 <- bedorfile1;
	bedfile2 <- bedorfile2;

	if(class(bedorfile1)=="data.frame")
	{
		bedfile1 = tempfile();
		write.table(bedorfile1, file=bedfile1, quote=F, col.names=F, row.names=F, sep="\t");
	}

	if(class(bedorfile2)=="data.frame")
	{
		bedfile2 = tempfile();
		write.table(bedorfile2, file=bedfile2, quote=F, col.names=F, row.names=F, sep="\t");
	}


	tbx <- NULL;
	#Only report those entries in A that have no overlap in B. Restricted by -f and -r.
	try(tbx <- read.table(pipe(paste("bedtools intersect -a ", bedfile1, "-b ", bedfile2 ))));

	if(class(bedorfile1)=="data.frame") unlink(bedfile1);
	if(class(bedorfile2)=="data.frame") unlink(bedfile2);

	return(tbx);
}

jaccard <-function( bedorfile1, bedorfile2 )
{
	bedfile1 <- bedorfile1;
	bedfile2 <- bedorfile2;

	if(class(bedorfile1)=="data.frame")
	{
		bedfile1 = tempfile();
		write.table(bedorfile1, file=bedfile1, quote=F, col.names=F, row.names=F, sep="\t");
	}

	if(class(bedorfile2)=="data.frame")
	{
		bedfile2 = tempfile();
		write.table(bedorfile2, file=bedfile2, quote=F, col.names=F, row.names=F, sep="\t");
	}

	tbx <- NULL;
	#Only report those entries in A that have no overlap in B. Restricted by -f and -r.
#system(paste("head ", bedfile2) )
	try(tbx <- read.table(pipe(paste("bedtools jaccard  -a ", bedfile1, "-b ", bedfile2 )), head=T));

	if(class(bedorfile1)=="data.frame") unlink(bedfile1);
	if(class(bedorfile2)=="data.frame") unlink(bedfile2);

	return(tbx);
}


dhs_overlap_removek562_H3k27ac <- function(peak_G1)
{
	tb.peak <- read.table(peak_G1);

	dhs_folder <- "/workdir/tc532/manuscript_code/database/merged/min150_sorted_data"
	k562_dnse_site <- "K562.bed";
	tb0 <- unintersect( peak_G1, paste(dhs_folder, k562_dnse_site, sep="/" ));
	tb0 <- unintersect( tb0, file.H3k27ac );

	df <- data.frame(dhs=list.files(dhs_folder, pattern="*.bed"), count=0, tb1=0, tb2=0, overlap=NA, ratio=0, jaccard=NA );

	for(i in 1:NROW(df))
	{
cat(paste(dhs_folder, df[i,1], sep="/" ), "\n");
		tb1 <- unintersect(  paste(dhs_folder, df[i,1], sep="/" ), paste(dhs_folder, k562_dnse_site, sep="/" ) );
		tb1 <- unintersect( tb1, file.H3k27ac );

		df[i,2] <- NROW(read.table(paste(dhs_folder, df[i,1], sep="/" )));
		df[i,3] <- NROW(tb1);
		df[i,4] <- NROW(tb0);

		if (NROW(tb0)==0 || NROW(tb1)==0 )
			next;

		tbx <- intersect(tb0, tb1);
		df[i,5] <- NROW(tbx);
		df[i,6] <- NROW(tbx)/NROW(tb.peak)
		tbj <- jaccard(  tb0, tb1 );
		if(NROW(tbj)>=1)
			df[i,7] <- tbj[1,"jaccard"];
	}

	return( df[ df$dhs != "K562.bed", ] );
}

if(1)
{
	df.G1 <- dhs_overlap_removek562_H3k27ac("../new-rf-201803/G1/G1.dREG.peak.score.bed.gz");
	df.G2 <- dhs_overlap_removek562_H3k27ac("../new-rf-201803/G2/G2.dREG.peak.score.bed.gz");
	df.G3 <- dhs_overlap_removek562_H3k27ac("../new-rf-201803/G3/G3.dREG.peak.score.bed.gz");
	df.G4 <- dhs_overlap_removek562_H3k27ac("../new-rf-201803/G4/G4.dREG.peak.score.bed.gz");
	df.G5 <- dhs_overlap_removek562_H3k27ac("../new-rf-201803/G5/G5.dREG.peak.score.bed.gz");
	df.G6 <- dhs_overlap_removek562_H3k27ac("../new-rf-201803/G6/G6.dREG.peak.score.bed.gz");
	df.G7 <- dhs_overlap_removek562_H3k27ac("../new-rf-201803/G7/G7.dREG.peak.score.bed.gz");

	save(df.G1, df.G2, df.G3, df.G4, df.G5, df.G6, df.G7, file="novel-dreg-removek562-h3k27ac.rdata");
}


dhs_overlap_removek562 <- function(peak_G1)
{
	tb.peak <- read.table(peak_G1);

	k562_dnse_site <- "K562.bed";
	tb0 <- unintersect( peak_G1, paste(dhs_folder, k562_dnse_site, sep="/" ));

	df <- data.frame(dhs=list.files(dhs_folder, pattern="*.bed"), count=0, tb1=0, tb2=0, overlap=NA, ratio=0, jaccard=NA );

	for(i in 1:NROW(df))
	{
		tb1 <- unintersect(  paste(dhs_folder, df[i,1], sep="/" ), paste(dhs_folder, k562_dnse_site, sep="/" ) );
		tbx <- intersect(tb0, tb1);
		df[i,2] <- NROW(read.table(paste(dhs_folder, df[i,1], sep="/" )));
		df[i,3] <- NROW(tb1);
		df[i,4] <- NROW(tb0);
		df[i,5] <- NROW(tbx);
		df[i,6] <- NROW(tbx)/NROW(tb.peak)
		tbj <- jaccard(  tb0, tb1 );
		if(NROW(tbj)>=1)
			df[i,7] <- tbj[1,"jaccard"];
	}

	return( df[ df$dhs != "K562.bed", ] );
}

if(1)
{
	df.G1 <- dhs_overlap_removek562("../new-rf-201803/G1/G1.dREG.peak.score.bed.gz");
	df.G2 <- dhs_overlap_removek562("../new-rf-201803/G2/G2.dREG.peak.score.bed.gz");
	df.G3 <- dhs_overlap_removek562("../new-rf-201803/G3/G3.dREG.peak.score.bed.gz");
	df.G4 <- dhs_overlap_removek562("../new-rf-201803/G4/G4.dREG.peak.score.bed.gz");
	df.G5 <- dhs_overlap_removek562("../new-rf-201803/G5/G5.dREG.peak.score.bed.gz");
	df.G6 <- dhs_overlap_removek562("../new-rf-201803/G6/G6.dREG.peak.score.bed.gz");
	df.G7 <- dhs_overlap_removek562("../new-rf-201803/G7/G7.dREG.peak.score.bed.gz");

	save(df.G1, df.G2, df.G3, df.G4, df.G5, df.G6, df.G7, file="novel-dreg-removek562.rdata");
}

dhs_overlap_keepk562 <- function(peak_G1)
{
	tb.peak <- read.table(peak_G1);

	df <- data.frame(dhs=list.files(dhs_folder, pattern="*.bed"), intersection=0, union=0, n_intersections=0, jaccard=NA );
	for(i in 1:NROW(df))
	{
		tbj <- jaccard(  paste(dhs_folder, df[i,1], sep="/" ), tb.peak );

		if(NROW(tbj)>=1)
		{
			df[i,2] <- tbj[1,"intersection"]
			df[i,3] <- tbj[1,"union.intersection"]
			df[i,4] <- tbj[1,"n_intersections"]
			df[i,5] <- tbj[1,"jaccard"]
		}
	}

	return( df);
}

if(1)
{
	df.G1 <- dhs_overlap_keepk562("../new-rf-201803/G1/G1.dREG.peak.score.bed.gz");
	df.G2 <- dhs_overlap_keepk562("../new-rf-201803/G2/G2.dREG.peak.score.bed.gz");
	df.G3 <- dhs_overlap_keepk562("../new-rf-201803/G3/G3.dREG.peak.score.bed.gz");
	df.G4 <- dhs_overlap_keepk562("../new-rf-201803/G4/G4.dREG.peak.score.bed.gz");
	df.G5 <- dhs_overlap_keepk562("../new-rf-201803/G5/G5.dREG.peak.score.bed.gz");
	df.G6 <- dhs_overlap_keepk562("../new-rf-201803/G6/G6.dREG.peak.score.bed.gz");
	df.G7 <- dhs_overlap_keepk562("../new-rf-201803/G7/G7.dREG.peak.score.bed.gz");

	save(df.G1, df.G2, df.G3, df.G4, df.G5, df.G6, df.G7, file="novel-dreg-keep562.rdata");
}



