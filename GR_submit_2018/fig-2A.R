make_big_bed<-function(file.nonneg )
{
	system( paste("cat ../k562/wgEncodeOpenChromDnaseK562PkV2.narrowPeak ../k562/GSM646567_hg19_wgEncodeUwDgfK562Pk.narrowPeak.txt ../k562/GSM646567_hg19_wgEncodeUwDgfK562Pk.macs2.narrowPeak | awk -v OFS='\\t' '{print $1,$2,$3}' - > ", file.nonneg ), intern=TRUE );

	nonneg_bed <- unique(read.table( pipe(paste("cat ", file.nonneg, " | sort-bed ", file.nonneg," | bedtools merge -i - -d 100 ", sep=" " ) ) )[,c(1:3)]);
	#nonneg_bed[,2] <- nonneg_bed[,2] - 100
	#idx.mis <- which(nonneg_bed[,2]<0);
	#if(length(idx.mis)>0) nonneg_bed[idx.mis,2] <- 0;
	#nonneg_bed[,3] <- nonneg_bed[,3] + 100
	write.table( nonneg_bed, file=file.nonneg, quote=F, row.name=F, col.names=F, sep="\t" );
}


file.K562.DHS.peak  <- "./k562.merged.broad.peak.bed";
make_big_bed( file.K562.DHS.peak );

tb.dREG <- read.table("/workdir/zw355/proj/prj10-dreg/new-rf-201803/G1/G1.dREG.peak.full.bed.gz");
tb.dNase <- read.table(file.K562.DHS.peak);
tb.H3K27ac <- read.table("/fs/cbsudanko/storage/data/hg19/k562/histones/wgEncodeBroadHistoneK562H3k27acStdAln.bed.gz");

cat("mean dREG=", mean(tb.dREG[,3]-tb.dREG[,2]), "dNase=", mean(tb.dNase[,3]-tb.dNase[,2]),"H3K27ac=", mean(tb.H3K27ac[,3]-tb.H3K27ac[,2]) , "\n");
cat("median dREG=", median(tb.dREG[,3]-tb.dREG[,2]), "dNase=", median(tb.dNase[,3]-tb.dNase[,2]),"H3K27ac=", median(tb.H3K27ac[,3]-tb.H3K27ac[,2]) , "\n");

df <- rbind(data.frame(src="dNase", pwidth=tb.dNase[,3]-tb.dNase[,2]),
      data.frame(src="dREG", pwidth=tb.dREG[,3]-tb.dREG[,2]),
      data.frame(src="H3K27ac", pwidth=tb.H3K27ac[,3]-tb.H3K27ac[,2]) )


library(plyr)
mu <- ddply(df, "src", summarise, grp.mean=mean(pwidth))
head(mu)

library(ggplot2);
library(plotly);

p <- ggplot(df, aes(pwidth, fill = src)) + geom_density(alpha = 0.2, size=0.15) +
        ylim( 0, 0.01 )+
        xlim( 0, 8000 )+
        xlab("Peak width") +
        ylab("Density") +
        scale_fill_discrete(name = "Data set") + 
        theme(axis.text.y=element_text(angle=90, hjust=1))+
        theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), 
               axis.line = element_line(colour = "black", size=0.3), legend.position = c(0.8, 0.78) )

p <- ggplotly(p)
ggsave("fig-2A.pdf", device="pdf", width=4, height=4)
