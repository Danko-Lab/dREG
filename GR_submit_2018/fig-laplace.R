
file.G1.dreg <- "../new-rf-201803/G1/G1.dREG.infp.bed.gz"
file.G1.dreg.peak <- "../new-rf-201803/G1/G1.dREG.peak.score.bed.gz"
file.G1.old<-"../old-model/test-G1/out.dREG.pred.gz"
file.G1.old.peak<-"../old-model/test-G1/out.dREG.peak.gz"

library(rmutil)
source("../dregHD-train/dREG_prob.R");

draw_hist<-function(file.pred, title, breaks=250, xlim=c(-0.2, 0.2), ylim=c(0, 50))
{
	if(class(file.pred)!="data.frame")
		tb.pred <- read.table( file.pred, header=F)
	else
		tb.pred <- file.pred;

	ypred <- tb.pred[,4]
	y.minus <- ypred [ypred<0];
	sigma.mn <- get_eps_laplace_sigma( c( y.minus, abs(y.minus), ypred[ypred==0]))
cat("sigma.mn=", sigma.mn, "\n");
	hist(ypred, breaks=breaks, prob=TRUE, main=title, xlab="dREG score", xlim=xlim, ylim=ylim, border="gray", cex.axis=1.5, cex.lab=1.5  );

	yfit <- dlaplace(seq(xlim[1], xlim[2], 0.01/2), m=0, s=sigma.mn, log=FALSE)
	lines(seq(xlim[1], xlim[2], 0.01/2), yfit, col = "blue", lwd = 1.2 )

	#remove normal distribution
	#yfit <- dnorm(seq(-1, 1, 0.02), mean=0, sd=sd(ypred), log=FALSE)
	#lines(seq(-1, 1, 0.02), yfit, col = "red", lwd = 1)

}

options(scipen =99) # not to use scientific notation when writing out

load("../dreg-train/k562.negative.bed.rdata");

#only located in Negative regions
file.neg = tempfile(fileext=".bed");
write.table(negative_bed, file.neg, quote=F, row.names=F, col.names=F, sep="\t")
tb.neg.new <- read.table(pipe(paste("zcat ", file.G1.dreg, " | bedtools intersect -b ", file.neg, "-a - -wb")))

file.neg.new = tempfile(fileext=".bed");
write.table(tb.neg.new, file=file.neg.new, quote=F, row.names=F, col.names=F, sep="\t")

#remove the peaks in Negative regions but detected by dREG
tb<-read.table(file.G1.dreg.peak);
tb[,2] <- tb[,2]-1000
tb[ tb[,2]<0, 2 ] <- 0;
tb[,3] <- tb[,3]+1000
file.temp.dreg = tempfile(fileext=".bed");
write.table(tb, file.temp.dreg, quote=F, row.names=F, col.names=F, sep="\t")

file.temp.dreg2 = tempfile(fileext=".bed");
system(paste("bedtools merge -i ", file.temp.dreg, "| sort-bed - >", file.temp.dreg2))
tb.pos.new <- read.table(pipe(paste("cat ", file.neg.new, " | bedtools intersect -b ", file.temp.dreg2, "-a - -wb")))
tb.pos.new <- unique( tb.pos.new[,c(1:4)] )

idx.pos <- match(paste(tb.pos.new[,1], tb.pos.new[,2], sep=":"), paste(tb.neg.new[,1], tb.neg.new[,2], sep=":"))
idx.pos <- idx.pos[!is.na(idx.pos)]


tb.new2 <- tb.neg.new[-idx.pos,][,c(1:5)];
#tb.new <- rbind( tb.new2[tb.new2[,4]>0, ], tb.neg.new[tb.neg.new[,4]<=0, c(1:5)] )
tb.new <- tb.new2


tb<-read.table(file.G1.dreg);
tb.neg <- tb[ which(tb[,4]<=0),]
ypred <- c(tb.neg[,4], abs(tb.neg[,4]));
sigma.mn <- get_eps_laplace_sigma( ypred );
cat("sigma.mn=", sigma.mn, "\n");

pdf("fig-laplace-5.pdf")

par(mar=c(5, 5, 4, 2) + 0.1)
draw_hist( tb.new[tb.new$V5==1,], breaks=400, xlim=c(-0.1, 0.1), ylim=c(0,60), " ");

xlim = c(-0.2, 0.2)
ylim = c(0, 50)
yfit <- dlaplace(seq(xlim[1], xlim[2], 0.01/2), m=0, s=sigma.mn, log=FALSE);
lines(seq(xlim[1], xlim[2], 0.01/2), yfit, col = "red", lwd = 2, lty=11);
legend("topright", c("All negative sites", "All negative values"), lty=c("solid", "11"), lwd=c(1,2), col=c("blue", "red"), cex=1.3);
dev.off();
