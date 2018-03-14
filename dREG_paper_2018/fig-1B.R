tpr.grocap.k562 <- c( 0.813, 0.414)
tpr.grocap.gm <- c( 0.961, 0.352)

tpr.h3kac.k562 <- c(0.584, 0.291)
tpr.h3kac.gm <- c(0.718, 0.225)

cols=c("#cb5b42", "#6985cd");

pdf("fig-1B.pdf");
x <- cbind(dREG1=tpr.grocap.k562, TFit1=tpr.grocap.gm, dREG2=tpr.h3kac.k562, TFit2=tpr.h3kac.gm);
colnames(x) <- c("K562/GRO-cap", "GM12878/GRO-cap", "K562/H3k27ac", "GM12878/H3k27ac");

barplot(x, xlab="", ylab= "TPR", beside=TRUE, col=cols, ylim=c(0,1), names.arg=c("K562/GRO-cap", "GM12878/GRO-cap", "K562/H3k27ac", "GM12878/H3k27ac"), cex.axis=0.95, cex.names=0.95, cex.lab=0.95)

# Place the legend at the top-left corner with no frame
# using rainbow colors
legend("topleft", c("dREG","TFit"), cex=1.2,  bty="n", fill=cols);

dev.off()