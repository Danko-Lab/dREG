
##
ss <- read.table("superset.nums.bed.gz")

## Using these values to construct a Venn diagriam by hand.
NROW(ss) # 160755

## dREG is ss[,4]
sum(ss[,4] == 1 & ss[,5] == 1 & ss[,6] == 1) # 23594
sum(ss[,4] == 1 & ss[,5] == 0 & ss[,6] == 0) # 4230
sum(ss[,4] == 1 & ss[,5] == 1 & ss[,6] == 0) # 394
sum(ss[,4] == 1 & ss[,5] == 0 & ss[,6] == 1) # 4148

## DNASE-1 is ss[,5]
#sum(ss[,4] == 1 & ss[,5] == 1 & ss[,6] == 1) ## REPEAT 
sum(ss[,4] == 0 & ss[,5] == 1 & ss[,6] == 0) # 9304
#sum(ss[,4] == 1 & ss[,5] == 1 & ss[,6] == 0) ## REPEAT.
sum(ss[,4] == 0 & ss[,5] == 1 & ss[,6] == 1) # 32129

## chromHMM is ss[,6]
#sum(ss[,4] == 1 & ss[,5] == 1 & ss[,6] == 1) ## REPEAT.
sum(ss[,4] == 0 & ss[,5] == 0 & ss[,6] == 1) # 86956
#sum(ss[,4] == 0 & ss[,5] == 1 & ss[,6] == 1) ## REPEAT.
#sum(ss[,4] == 1 & ss[,5] == 0 & ss[,6] == 1) ## REPEAT.

## Insulator is ss[,7]
sum(ss[,4] == 1 & ss[,5] == 1 & ss[,6] == 1 & ss[,7] == 1) # 1838
sum(ss[,4] == 0 & ss[,5] == 0 & ss[,6] == 1 & ss[,7] == 1) # 22481
sum(ss[,4] == 0 & ss[,5] == 1 & ss[,6] == 1 & ss[,7] == 1) # 16291
sum(ss[,4] == 1 & ss[,5] == 0 & ss[,6] == 1 & ss[,7] == 1) # 216

###############################################
## H3K27ac is ss[,8]
sum(ss[,4] == 1 & ss[,8] == 1 & ss[,5] == 1) # 17385
sum(ss[,4] == 0 & ss[,8] == 1 & ss[,5] == 1) # 3541
sum(ss[,4] == 1 & ss[,8] == 0 & ss[,5] == 1) # 5769
5769/(3541+5769) ##62% more likely to have pol ii than h3k27ac
## A larger fraction have JUST dREG (as opposed to JUST H3K27ac).
## Does this support a model of transcription prior to H3K27 acetylation?
## Alternatively, it might be a sensitivity issue with H3K27ac, or a high FPR for dREG.
## Requiring DNAse-1 as a tie breaker equalizes these... So this reflects dREG FPR?!
## Anyway, nothing conclusive.

###############################################
## Express as a fraction.

## dREG is ss[,4]
sum(ss[,4] == 1 & ss[,5] == 1 & ss[,6] == 1)/NROW(ss) # 0.1467699
sum(ss[,4] == 1 & ss[,5] == 0 & ss[,6] == 0)/NROW(ss) # 0.02631333
sum(ss[,4] == 1 & ss[,5] == 1 & ss[,6] == 0)/NROW(ss) # 0.002450935
sum(ss[,4] == 1 & ss[,5] == 0 & ss[,6] == 1)/NROW(ss) # 0.02580324

## DNASE-1 is ss[,5]
#sum(ss[,4] == 1 & ss[,5] == 1 & ss[,6] == 1)/NROW(ss) ## REPEAT 
sum(ss[,4] == 0 & ss[,5] == 1 & ss[,6] == 0)/NROW(ss) # 0.05787689
#sum(ss[,4] == 1 & ss[,5] == 1 & ss[,6] == 0)/NROW(ss) ## REPEAT.
sum(ss[,4] == 0 & ss[,5] == 1 & ss[,6] == 1)/NROW(ss) # 0.1998631

## chromHMM is ss[,6]
#sum(ss[,4] == 1 & ss[,5] == 1 & ss[,6] == 1)/NROW(ss) ## REPEAT.
sum(ss[,4] == 0 & ss[,5] == 0 & ss[,6] == 1)/NROW(ss) # 0.5409225
#sum(ss[,4] == 0 & ss[,5] == 1 & ss[,6] == 1)/NROW(ss) ## REPEAT.
#sum(ss[,4] == 1 & ss[,5] == 0 & ss[,6] == 1)/NROW(ss) ## REPEAT.

## Insulator is ss[,7]
sum(ss[,4] == 0 & ss[,5] == 0 & ss[,6] == 1 & ss[,7] == 1)/NROW(ss) # 0.1398464
sum(ss[,4] == 0 & ss[,5] == 1 & ss[,6] == 1 & ss[,7] == 1)/NROW(ss) # 0.1013405
sum(ss[,4] == 1 & ss[,5] == 0 & ss[,6] == 1 & ss[,7] == 1)/NROW(ss) # 0.00134366

sum(ss[,4] == 1 & ss[,5] == 1 & ss[,6] == 0 & ss[,7] == 0)/NROW(ss) # 0.002450935

#For Fig. 1 Venn diagriam
sum(ss[,4] == 0 & ss[,5] == 0 & ss[,6] == 1 & ss[,7] == 0)/NROW(ss) # 0.4010762
sum(ss[,4] == 0 & ss[,6] == 1 & ss[,7] == 1)/NROW(ss) # 0.2411869
sum(ss[,4] == 0 & ss[,5] == 1 & ss[,6] == 1 & ss[,7] == 0)/NROW(ss) # 0.0985226
sum(ss[,4] == 1 & ss[,5] == 1 & ss[,6] == 1 & ss[,7] == 0)/NROW(ss) # 0.1353364

sum(ss[,4] == 0 & ss[,5] == 1 & ss[,6] == 0 & ss[,7] == 0)/NROW(ss) # 0.05787689
sum(ss[,4] == 1 & ss[,5] == 0 & ss[,6] == 0 & ss[,7] == 0)/NROW(ss) # 0.02631333
sum(ss[,4] == 1 & ss[,5] == 0 & ss[,6] == 1 & ss[,7] == 0)/NROW(ss) # 0.02445958
sum(ss[,4] == 1 & ss[,5] == 1 & ss[,6] == 1 & ss[,7] == 1)/NROW(ss) # 0.01143355


#################################################################
## ADD GM12878 AND:
## Use a stacked barplot?!  Allows us to show both k562 and gm12878.
ss <- read.table("superset.nums.bed")
sg <- read.table("gm12878.superset.nums.bed")
sc <- read.table("cd4.superset.nums.bed")
sh <- read.table("hela.superset.nums.bed")

ka <- c("chromHMM Only"= sum(ss[,4] == 0 & ss[,5] == 0 & ss[,6] == 1 & ss[,7] == 0)/ NROW(ss),
  "chromHMM & DNAse-1"= sum(ss[,4] == 0 & ss[,5] == 1 & ss[,6] == 1 & ss[,7] == 0)/ NROW(ss),
  "chromHMM & DNAse-1 & dREG"= sum(ss[,4] == 1 & ss[,5] == 1 & ss[,6] == 1)/ NROW(ss), 
  "Insulator"= sum(ss[,7] == 1 & ss[,4] == 0)/ NROW(ss),
  "Other"= sum((ss[,4] == 0 & ss[,5] == 1 & ss[,6] == 0) | (ss[,4] == 1 & ss[,5] == 0 & ss[,6] == 1) | 
			(ss[,4] == 1 & ss[,5] == 0 & ss[,6] == 0) | (ss[,4] == 1 & ss[,5] == 1 & ss[,6] == 0) )/ NROW(ss))
			
ga <- c("chromHMM Only"= sum(sg[,4] == 0 & sg[,5] == 0 & sg[,6] == 1 & sg[,7] == 0)/ NROW(sg),
  "chromHMM & DNAse-1"= sum(sg[,4] == 0 & sg[,5] == 1 & sg[,6] == 1 & sg[,7] == 0)/ NROW(sg),
  "chromHMM & DNAse-1 & dREG"= sum(sg[,4] == 1 & sg[,5] == 1 & sg[,6] == 1)/ NROW(sg), 
  "Insulator"= sum(sg[,7] == 1 & sg[,4] == 0)/ NROW(sg),
  "Other"= sum((sg[,4] == 0 & sg[,5] == 1 & sg[,6] == 0) | (sg[,4] == 1 & sg[,5] == 0 & sg[,6] == 1) | 
			(sg[,4] == 1 & sg[,5] == 0 & sg[,6] == 0) | (sg[,4] == 1 & sg[,5] == 1 & sg[,6] == 0) )/ NROW(sg))


ca <- c("chromHMM Only"= sum(sc[,4] == 0 & sc[,5] == 0 & sc[,6] == 1 & sc[,7] == 0)/ NROW(sc),
  "chromHMM & DNAse-1"= sum(sc[,4] == 0 & sc[,5] == 1 & sc[,6] == 1 & sc[,7] == 0)/ NROW(sc),
  "chromHMM & DNAse-1 & dREG"= sum(sc[,4] == 1 & sc[,5] == 1 & sc[,6] == 1)/ NROW(sc), 
  "Insulator"= sum(sc[,7] == 1 & sc[,4] == 0)/ NROW(sc),
  "Other"= sum((sc[,4] == 0 & sc[,5] == 1 & sc[,6] == 0) | (sc[,4] == 1 & sc[,5] == 0 & sc[,6] == 1) | 
			(sc[,4] == 1 & sc[,5] == 0 & sc[,6] == 0) | (sc[,4] == 1 & sc[,5] == 1 & sc[,6] == 0) )/ NROW(sc))
			

ha <- c("chromHMM Only"= sum(sh[,4] == 0 & sh[,5] == 0 & sh[,6] == 1 & sh[,7] == 0)/ NROW(sh),
  "chromHMM & DNAse-1"= sum(sh[,4] == 0 & sh[,5] == 1 & sh[,6] == 1 & sh[,7] == 0)/ NROW(sh),
  "chromHMM & DNAse-1 & dREG"= sum(sh[,4] == 1 & sh[,5] == 1 & sh[,6] == 1)/ NROW(sh), 
  "Insulator"= sum(sh[,7] == 1 & sh[,4] == 0)/ NROW(sh),
  "Other"= sum((sh[,4] == 0 & sh[,5] == 1 & sh[,6] == 0) | (sh[,4] == 1 & sh[,5] == 0 & sh[,6] == 1) | 
			(sh[,4] == 1 & sh[,5] == 0 & sh[,6] == 0) | (sh[,4] == 1 & sh[,5] == 1 & sh[,6] == 0) )/ NROW(sh))
						
			
df <- data.frame(k562= ka, gm12878= ga, cd4= ca, hela= ha)
df

require(ggplot2)

qplot(sb, data=df, geom="bar", fill=factor(label))

####################################################
## And also check CD4+ T-cells (for good measure).
a <- read.table("cd4.superset.nums.bed")
sum(a$V4 == 1 & a$V5 == 1 & a$V6 == 1)/NROW(a)
sum(a$V4 == 0 & a$V5 == 1 & a$V6 == 1)/NROW(a)
sum(a$V4 == 0 & a$V5 == 0 & a$V6 == 1)/NROW(a)
sum(a$V4 == 0 & a$V5 == 1 & a$V6 == 0)/NROW(a)
######

######################################################
## Print overlaps for Fig. 3B.
library(lattice)

x_cex=y_cex=1.5
x_lab_cex=y_lab_cex=1.55
col=c("dark green", "dark blue")

overlap <- read.csv(text="Mark, Overlap
DNAse-1, 27
chromHMM, 22
H3K27ac, 72
H3K9ac, 86", header=TRUE)

barchart(Overlap~Mark,data=overlap, col="dark blue",
         scales=list(x=list(rot=65,cex=x_cex), y=list(rot=0, cex=y_cex)), ylim=c(0,100), ylab="Fraction of Sites Transcribed",
		 auto.key=list(columns = 2, space = "top"),
		 lwd=1, pch=1, xlab= list(cex=x_lab_cex))
