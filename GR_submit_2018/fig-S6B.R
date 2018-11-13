tpr.grocap.G8 <- c( 0.819, 0.821, 0.802)
tpr.grocap.gm <- c( 0.961, 0.966, 0.966)
tpr.h3k27ac.G8 <- c( 0.645, 0.626, 0.605)
tpr.h3k27ac.gm <- c( 0.718, 0.719, 0.702)
tpr.h3k27ac.gh <- c( 0.841, 0.804, 0.803)


cols=c("#cb5b42", "#6985cd", "#cbf21f");

pdf("fig-S6B.pdf");
x <- cbind( C11=tpr.grocap.G8,  C12=tpr.grocap.gm,  #C14=tpr.grocap.gh, 
            C21=tpr.h3k27ac.G8, C22=tpr.h3k27ac.gm, C23=tpr.h3k27ac.gh);
x.names <- c("G8/GRO-cap", "GM/GRO-cap", "G8/H3K27ac", "GM/H3K27ac", "GH/H3K27ac");

barplot(x, xlab="", ylab= "TPR", beside=TRUE, col=cols, ylim=c(0,1.05), names.arg=x.names, cex.axis=0.95, cex.names=0.8, cex.lab=0.8)

# Place the legend at the top-left corner with no frame
# using rainbow colors
legend("topright", c("Current Model(Trained G1,2,3,5,6)","Trained by G1,2,3,5,8", "Trained by G1,2,3,5,M"), cex=1.0,  bty="n", fill=cols);

dev.off()