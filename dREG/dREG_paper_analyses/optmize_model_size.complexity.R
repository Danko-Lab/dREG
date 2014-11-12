## For individual TSS.
system("sort-bed /usr/projects/GROseq.parser/tss_new/hg19.k562.new_hmm2.bed > ~/tmp.bed")

## Combine into clusters.
system("cat /usr/projects/GROseq.parser/tss_new/hg19.k562.new_hmm2.bed | grep -v \"^t\" | sort-bed - | awk 'BEGIN {OFS=\"\t\"} {print $1,$2,$3,\"N\",$5,\"+\"}' > andrehmm.nostrand.bed")
system("~/bin/tcolapse.tss andrehmm.nostrand.bed 0 1000 1 | sort-bed - > andrehmm.nostrand.merge.bed")

#system("cat predictions.bed | awk 'BEGIN{OFS=\"\t\"} {print $1,$2,$3,\"N\",$4}' | bedmap --header --max ~/tmp.bed - | grep \"NAN\" -c", intern=TRUE)

 require(featureDetector)
 setwd("/usr/projects/GROseq.parser/tss_detecter/")

 ps_plus_path  <- "/usr/data/GROseq.parser/hg19/k562/proseq/K562_unt.subsamp10pct.bed.gz_plus.bw"
 ps_minus_path <- "/usr/data/GROseq.parser/hg19/k562/proseq/K562_unt.subsamp10pct.bed.gz_minus.bw"

##############################################################################
## Using both strands. + and - strand ...
makeUseAND <- function(window=100, depth=0, step=100) {
 inf_positions_both <- get_informative_positions(ps_plus_path, ps_minus_path, window= window, depth= depth, step=step, use_OR=FALSE) 
 final_data <- data.frame(inf_positions_both)
 options(scipen=25)
 write.table(final_data, file="tmp.both.values.bed", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
 NROW(final_data)

 n_missed <- system(paste("cat tmp.both.values.bed | awk 'BEGIN{OFS=\"\t\"} {print $1,$2-",ceiling(step/2),",$3+",ceiling(step/2),"}' | bedmap --count andrehmm.nostrand.merge.bed - | grep \"^0$\" -c", sep=""), intern=TRUE)
 total <- system("grep \"\" -c andrehmm.nostrand.merge.bed", intern=TRUE)
 print(paste(n_missed, total))
 frac <- as.double(n_missed)/as.double(total)
 
 return(list(NROW(final_data), frac))
}

##############################################################################
## If we look near every position w/ >threshold reads ...
makeUseOR <- function(window=100, depth=1, step=100) {

 inf_positions <- get_informative_positions(ps_plus_path, ps_minus_path, window= window, depth= depth, step=step, use_OR=TRUE) 
 final_data <- data.frame(inf_positions)
 options(scipen=25)
 write.table(final_data, file="tmp.values.bed", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
 NROW(final_data)

 ## Conceptually, if we looked over this range...
 n_missed <- system(paste("cat tmp.values.bed | awk 'BEGIN{OFS=\"\t\"} {print $1,$2-",ceiling(step/2),",$3+",ceiling(step/2),"}' | bedmap --count andrehmm.nostrand.merge.bed - | grep \"^0$\" -c", sep=""), intern=TRUE)
 total <- system("grep \"\" -c andrehmm.nostrand.merge.bed", intern=TRUE)
 print(paste(n_missed, total))
 frac <- as.double(n_missed)/as.double(total)
 
 return(list(NROW(final_data), frac))
}

##############################################################################
## Combined ... might be a good tradeoff ...
makeUseANDOR <- function(window=500, depth=0, step=100) {
 inf_positions_both <- get_informative_positions(ps_plus_path, ps_minus_path, window= window, depth= depth, step=step, use_ANDOR=TRUE) 
 final_data <- data.frame(inf_positions_both)
 options(scipen=25)
 write.table(final_data, file="tmp.ANDOR.values.bed", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
 NROW(final_data)

 n_missed <- system(paste("cat tmp.ANDOR.values.bed | awk 'BEGIN{OFS=\"\t\"} {print $1,$2-",ceiling(step/2),",$3+",ceiling(step/2),"}' | bedmap --count andrehmm.nostrand.merge.bed - | grep \"^0$\" -c", sep=""), intern=TRUE)
 total <- system("grep \"\" -c andrehmm.nostrand.merge.bed", intern=TRUE)
 print(paste(n_missed, total))
 frac <- as.double(n_missed)/as.double(total)
 
 return(list(NROW(final_data), frac))
}

plot_data <- function(nreads, sens, med, x_lab) {
options(scipen=4)

par(mar=c(5, 4, 2, 8) + 0.1)

plot(nreads, sens, axes=F, ylim=c(0,max(sens)), xlab="", ylab="",type="b",col="black", main="")
points(nreads,sens,pch=20,col="black")
axis(2, ylim=c(0,max(sens)),col="black",lwd=2)
mtext(2,text="Number of missed TSS.",line=2)

axis(1, ylim=c(0,max(nreads)),col="black",lwd=2)
mtext(1,text=x_lab,line=2)

par(new=T)

plot(nreads, med, axes=F, ylim=c(0,max(med)), col="dark red", xlab="", ylab="", type="b",lty=2, main="",lwd=2)
axis(4, ylim=c(0,max(med)),lwd=2,line=3.5, col="dark red")
points(nreads, med,pch=20, col="dark red")
mtext(4,text="Threshold number of reads.",line=5.5, col="dark red")

}

pdf("TestOptDepthPlots.pdf")

## Test each systematically ...
n_miss <- integer(0)
compl  <- integer(0)
rdepth <- c(0:4)
for(depth in rdepth) {
 a <- makeUseAND(window=1000, depth=depth, step=100)
 n_miss <- c(n_miss, a[[2]])
 compl  <- c(compl,  as.integer(a[[1]]))
}
plot_data(rdepth, n_miss, compl, "Number of positions to evaluate.")

n_miss <- integer(0)
compl  <- integer(0)
steps <- c(5,10,50,75,100,150,200,500)
for(step in steps) {
 a <- makeUseAND(window=1000, depth=0, step=step)
 n_miss <- c(n_miss, a[[2]])
 compl  <- c(compl,  as.integer(a[[1]]))
}
plot_data(steps, n_miss, compl, "Step size (window 1000)")


n_miss <- integer(0)
compl  <- integer(0)
rdepth <- c(0:5)
for(depth in rdepth) {
 makeUseOR(window=50, depth=depth, step=1000)
 n_miss <- c(n_miss, a[[2]])
 compl  <- c(compl,  as.integer(a[[1]]))
}
plot_data(rdepth, n_miss, compl, "Number of positions to evaluate.")

dev.off()

## Test a good combined version ...

makeUseAND(window=1000, depth=0, step=100)
makeUseOR(window=100, depth=1, step=100)
system("cat tmp.both.values.bed | awk 'BEGIN{OFS=\"\t\"} {print $1,($2-51),($3+51)}' > tmp.bed")
system("cat tmp.values.bed | awk 'BEGIN{OFS=\"\t\"} {print $1,($2-51),($3+51)}' >> tmp.bed")
n_missed <- system("cat tmp.bed | sort-bed - | bedmap --count andrehmm.nostrand.merge.bed - | grep \"^0\" -c", intern=TRUE)
total <- system("grep \"\" -c andrehmm.nostrand.merge.bed", intern=TRUE)
print(as.double(n_missed)/as.double(total))

makeUseANDOR()

## NOTE: This is wrong, because there's likely to be some overlap if implemented in featureDetector
#n_sites <- system("cat tmp.bed | grep \"\" -c", intern=TRUE)
#print(n_sites)

system("cat ~/tmp.bed | grep \"\" -c ")
system("rm ~/tmp.bed")


