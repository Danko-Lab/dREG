
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

	#yfit <- dlaplace(seq(xlim[1], xlim[2], 0.01/2), m=0, s=sigma.mn, log=FALSE)
	#lines(seq(xlim[1], xlim[2], 0.01/2), yfit, col = "blue", lwd = 1.2 )

	#remove normal distribution
	#yfit <- dnorm(seq(-1, 1, 0.02), mean=0, sd=sd(ypred), log=FALSE)
	#lines(seq(-1, 1, 0.02), yfit, col = "red", lwd = 1)

}

if(0)
{
	pdf("fig-S4-laplace-1.pdf")
	par(mfrow=c(1,2))
	draw_hist(file.G1.dreg, "new dREG", xlim=c(-0.5, 1.5), ylim=c(0,20));
	draw_hist(file.G1.old,  "old dREG", xlim=c(-0.5, 1.5), ylim=c(0,20));
	dev.off()
}

if(0)
{
	load("../dreg-train/k562.negative.bed.rdata");
	options(scipen =99) # not to use scientific notation when writing out
	file.neg = tempfile(fileext=".bed");

	write.table(negative_bed, file.neg, quote=F, row.names=F, col.names=F, sep="\t")
	tb.neg.new <- read.table(pipe(paste("zcat ", file.G1.dreg, " | bedtools intersect -b ", file.neg, "-a - -wb")))
	tb.neg.old <- read.table(pipe(paste("zcat ", file.G1.old, " | bedtools intersect -b ", file.neg, "-a - -wb")))

	tb.neg.new <- unique( tb.neg.new[,c(1:4)] )
	tb.neg.old <- unique( tb.neg.old[,c(1:4)] )

	pdf("fig-S4-laplace-2.pdf")
	par(mfrow=c(1,2))
	draw_hist( tb.neg.new, "new dREG");
	draw_hist( tb.neg.old, "old dREG");
	dev.off();
}


if(0)
{
options(scipen =99) # not to use scientific notation when writing out

#only located in Negative regions
file.neg = tempfile(fileext=".bed");
write.table(negative_bed, file.neg, quote=F, row.names=F, col.names=F, sep="\t")
tb.neg.new <- read.table(pipe(paste("zcat ", file.G1.dreg, " | bedtools intersect -b ", file.neg, "-a - -wb")))
tb.neg.old <- read.table(pipe(paste("zcat ", file.G1.old, " | bedtools intersect -b ", file.neg, "-a - -wb")))

file.neg.new = tempfile(fileext=".bed");
write.table(tb.neg.new, file=file.neg.new, quote=F, row.names=F, col.names=F, sep="\t")
file.neg.old = tempfile(fileext=".bed");
write.table(tb.neg.old, file=file.neg.old, quote=F, row.names=F, col.names=F, sep="\t")

#remove the peaks in Negative regions but detected by dREG
tb<-read.table(file.G1.dreg.peak);
tb[,2] <- tb[,2]-20000
tb[ tb[,2]<0, 2 ] <- 0;
tb[,3] <- tb[,3]+20000
file.temp.dreg = tempfile(fileext=".bed");
write.table(tb, file.temp.dreg, quote=F, row.names=F, col.names=F, sep="\t")

tb<-read.table(file.G1.old.peak);
tb[,2] <- tb[,2]-20000
tb[ tb[,2]<0, 2 ] <- 0;
tb[,3] <- tb[,3]+20000
file.temp.old = tempfile(fileext=".bed");
write.table(tb, file.temp.old, quote=F, row.names=F, col.names=F, sep="\t")

file.temp.dreg2 = tempfile(fileext=".bed");
system(paste("bedtools merge -i ", file.temp.dreg, "| sort-bed - >", file.temp.dreg2))
tb.pos.new <- read.table(pipe(paste("cat ", file.neg.new, " | bedtools intersect -b ", file.temp.dreg2, "-a - -wb")))

file.temp.old2 = tempfile(fileext=".bed");
system(paste("bedtools merge -i ", file.temp.old, "| sort-bed - >", file.temp.old2))
tb.pos.old <- read.table(pipe(paste("cat ", file.neg.old, " | bedtools intersect -b ", file.temp.old2, "-a - -wb")))

tb.pos.new <- unique( tb.pos.new[,c(1:4)] )
tb.pos.old <- unique( tb.pos.old[,c(1:4)] )

idx.pos <- match(paste(tb.pos.new[,1], tb.pos.new[,2], sep=":"), paste(tb.neg.new[,1], tb.neg.new[,2], sep=":"))
idx.pos <- idx.pos[!is.na(idx.pos)]
tb.new <- tb.neg.new[-idx.pos,]

idx.pos <- match(paste(tb.pos.old[,1], tb.pos.old[,2], sep=":"), paste(tb.neg.old[,1], tb.neg.old[,2], sep=":"))
idx.pos <- idx.pos[!is.na(idx.pos)]
tb.old <- tb.neg.old[-idx.pos,]

pdf("fig-S4-laplace-2.pdf")
par(mfrow=c(1,2))
draw_hist( tb.new[tb.new$V5==1,], "new dREG");
draw_hist( tb.old, "old dREG");
dev.off();

}

if(0)
{
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

pdf("fig-S4-laplace-3.pdf")
draw_hist( tb.new[tb.new$V5==1,], "new dREG");
dev.off();

}

if(0)
{
options(scipen =99) # not to use scientific notation when writing out

tb<-read.table(file.G1.dreg);
tb.neg <- tb[ which(tb[,4]<=0),]
ypred <- c(tb.neg[,4], abs(tb.neg[,4]));
sigma.mn <- get_eps_laplace_sigma( ypred );
cat("sigma.mn=", sigma.mn, "\n");
xlim = c(-0.2, 0.2)
ylim = c(0, 50)
title="Negative only"

pdf("fig-S4-laplace-4.pdf")
hist(ypred, breaks=250, prob=TRUE, main=title, xlab="dREG score", xlim=xlim, ylim=ylim  );
rect(0.001,0,0.2, 50, border="white", col="white");
yfit <- dlaplace(seq(xlim[1], xlim[2], 0.01), m=0, s=sigma.mn, log=FALSE);
lines(seq(xlim[1], xlim[2], 0.01), yfit, col = "blue", lwd = 1);
dev.off();

}


if(1)
{
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

sigma.all <- get_eps_laplace_sigma( tb.new[tb.new$V5==1,4] )

pdf("fig-S4-laplace-5.pdf")

par(mar=c(5, 5, 4, 2) + 0.1)
draw_hist( tb.new[tb.new$V5==1,], breaks=400, xlim=c(-0.1, 0.1), ylim=c(0,60), " ");

xlim = c(-0.2, 0.2)
ylim = c(0, 50)

yfit <- dlaplace(seq(xlim[1], xlim[2], 0.01/2), m=0, s=sigma.all, log=FALSE)
lines(seq(xlim[1], xlim[2], 0.01/2), yfit, col = "blue", lwd = 1.2 )

yfit <- dlaplace(seq(xlim[1], xlim[2], 0.01/2), m=0, s=sigma.mn, log=FALSE);
lines(seq(xlim[1], xlim[2], 0.01/2), yfit, col = "red", lwd = 2, lty=11);
legend("topright", c("All negative sites", "All negative values"), lty=c("solid", "11"), lwd=c(1,2), col=c("blue", "red"), cex=1.3);
dev.off();

}

if(1)
{
	sigma.all <- get_eps_laplace_sigma( tb.new[tb.new$V5==1,4] )

	yp <- plaplace(tb.new[tb.new$V5==1,4], m=0, s=sigma.mn);
	png("qqplot.neg.neg.png")
	qqplot( -log10(runif(NROW(yp))), -log10(yp),xlim=c(0, 10), ylim=c(0, 10), xlab="Theoretical quantile", ylab="All Negative sites" );
	segments(0, 0, 10,10);
	dev.off()

	yp <- plaplace(tb[tb$V5==1,4]  , m=0, s=sigma.mn);
	png("qqplot.neg.all.png")
	qqplot(-log10(runif(NROW(yp))), -log10(yp), xlim=c(0, 10), ylim=c(0, 10), xlab="Theoretical quantile", ylab="All sites" );
	segments(0, 0, 10,10);
	dev.off()


	yp <- plaplace(tb.new[tb.new$V5==1,4], m=0, s=sigma.all);
	png("qqplot.all.neg.png")
	qqplot( -log10(runif(NROW(yp))), -log10(yp),xlim=c(0, 10), ylim=c(0, 10), xlab="Theoretical quantile", ylab="All Negative sites" );
	segments(0, 0, 10,10);
	dev.off()

	yp <- plaplace(tb[tb$V5==1,4], m=0, s=sigma.all);
	png("qqplot.all.all.png")
	qqplot( -log10(runif(NROW(yp))), -log10(yp),xlim=c(0, 10), ylim=c(0, 10), xlab="Theoretical quantile", ylab="All sites" );
	segments(0, 0, 10,10);
	dev.off()


	yp1 <- plaplace(tb.new[tb.new$V5==1,4], m=0, s=sigma.mn);
	yp2 <- plaplace(tb.new[tb.new$V5==1,4], m=0, s=sigma.all);
	png("qqplot.neg.png")
	qqplot( -log10(yp2), -log10(yp1),xlim=c(0, 10), ylim=c(0, 10), xlab="Theoretical quantile", ylab="All negative sites" );
	segments(0, 0, 10,10);
	dev.off()


	yp1 <- plaplace(tb[tb$V5==1,4], m=0, s=sigma.mn);
	yp2 <- plaplace(tb[tb$V5==1,4], m=0, s=sigma.all);
	png("qqplot.all.png")
	qqplot( -log10(yp2), -log10(yp1),xlim=c(0, 10), ylim=c(0, 10), xlab="Theoretical quantile", ylab="All sites" );
	segments(0, 0, 10,10);
	dev.off()


	#yp1 <- plaplace(tb[tb$V5==1,4], m=0, s=sigma.mn);
	#yp2 <- rlaplace(NROW(tb[tb$V5==1,4]), m=0, s=sigma.mn);
	#png("qqplot.all.ran.png")
	#qqplot( -log10(yp2), -log10(yp1),xlim=c(0, 10), ylim=c(0, 10), xlab="Theoretical quantile", ylab="All sites" );
	#segments(0, 0, 10,10);
	#dev.off()


}