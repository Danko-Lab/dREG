tpr.grocap.G7 <- c( 0.813, 0.414)
tpr.grocap.G8 <- c( 0.819, 0.381)
tpr.grocap.gm <- c( 0.961, 0.352)
tpr.grocap.gh <- c( 0.432, 0.214)

tpr.h3kac.G7 <- c(0.584, 0.291)
tpr.h3kac.G8 <- c(0.654, 0.277)
tpr.h3kac.gm <- c(0.718, 0.225)
tpr.h3kac.gh <- c(0.849, 0.373)

cols=c("#cb5b42", "#6985cd");

pdf("fig-1B.pdf");
x <- cbind(C11=tpr.grocap.G7, C12=tpr.grocap.G8, C13=tpr.grocap.gm, #C14=tpr.grocap.gh, 
           C21=tpr.h3kac.G7,  C22=tpr.h3kac.G8,  C23=tpr.h3kac.gm,  C24=tpr.h3kac.gh);
x.names <- c("G7/GRO-cap", "G8", "GM12878", #"HCT116", 
             "G7/H3K27ac", "G8", "GM12878", "HCT116");

barplot(x, xlab="", ylab= "TPR", beside=TRUE, col=cols, ylim=c(0,1), names.arg=x.names, cex.axis=0.95, cex.names=0.6, cex.lab=0.95)

# Place the legend at the top-left corner with no frame
# using rainbow colors
legend("topleft", c("dREG","TFit"), cex=1.2,  bty="n", fill=cols);

dev.off()