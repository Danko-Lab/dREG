library(lattice)
library(boot)
library(featureDetector)

x_cex=y_cex=1.5
x_lab_cex=y_lab_cex=1.55
col=c("dark green", "dark blue")

setwd("/usr/projects/GROseq.parser/tss_detecter/")

##########################################
## Get data from dREG data.
load("asvm.RData")
gdm <- genomic_data_model(window_sizes= c(10, 25, 50, 500, 5000), half_nWindows= c(10, 10, 30, 20, 20))

ps_plus_path <- "/usr/data/GROseq.parser/hg19/k562/proseq/K562_unt.sort.bed.gz_plus.bw"
ps_minus_path <- "/usr/data/GROseq.parser/hg19/k562/proseq/K562_unt.sort.bed.gz_minus.bw" 
allSites_bed <- read.table("/usr/projects/GROseq.parser/tss_detecter/k562.allSites.bed")

n_eval <- 30000
inf_positions <- get_informative_positions(ps_plus_path, ps_minus_path, depth= 0, step=50, use_ANDOR=TRUE, use_OR=FALSE)
tset <- get_test_set(inf_positions, GROcap_tss_bed, n_eval, allow= allSites_bed, enrich_negative_near_pos= 0, extra_enrich_frac= 0, avoid_dist= 100)
remove(inf_positions)

pdf("fig1.pdf")

##########################################
## Fig. 1A (Extra panel for brwoser shot?!)
##########################################

##########################################
## Fig. 1B (Data model)
  plot.gdm(gdm, positive_positions= tset[c(1:(NROW(tset)/2)),], bigwig_plus= ps_plus_path, bigwig_minus= ps_minus_path)
  


##########################################
## Fig. 1C (ROC plot)

pred_val<- eval_reg_svm(gdm, asvm, tset, ps_plus_path, ps_minus_path, batch_size= 10000)
dREG <- logreg.roc.calc(tset[,4],pred_val)

AUC <- roc.auc(dREG)
AUC
## Get data from alternative methods.
roc.plot(dREG, main=paste("AUC=", AUC))

##########################################
## Fig. 1D (Accuracy on data types)
acc <- read.csv(text="Cell type, Comparison, Sensitivity
K562, GRO-cap Pairs, 93
K562, GRO-cap Sites, 84
K562, CAGE TSS, 92
GM12878, GRO-cap Pairs, 97
GM12878, GRO-cap Sites, 90
GM12878, CAGE TSS, 90", header=TRUE)

barchart(Comparison~Sensitivity,data=acc,groups=Cell.type, #col=col,
         scales=list(x=list(rot=0,cex=x_cex), y=list(rot=35, cex=y_cex)), xlim=c(0,100),
		 auto.key=list(columns = 2, space = "top"),
		 lwd=1, pch=1, xlab= list(cex=x_lab_cex))

		 
acc <- read.csv(text="Cell type, Comparison, Sensitivity
K562, GRO-cap Gene body, 85
K562, GRO-cap Enhancers, 78
K562, GRO-cap Promoters, 94
GM12878, GRO-cap Gene body, 90
GM12878, GRO-cap Enhancers, 85
GM12878, GRO-cap Promoters, 95", header=TRUE)

barchart(Comparison~Sensitivity,data=acc,groups=Cell.type, #col=col,
         scales=list(x=list(rot=0,cex=x_cex), y=list(rot=35, cex=y_cex)), xlim=c(0,100),
		 auto.key=list(columns = 2, space = "top"),
		 lwd=1, pch=1, xlab= list(cex=x_lab_cex))

		 
acc <- read.csv(text="Cell type, Method, Ernst classes, Sensitivity
K562, GRO-cap, 1_Active_Promoter, 80
K562, dREG, 1_Active_Promoter, 88
K562, GRO-cap, 2_Weak_Promoter, 22
K562, dREG, 2_Weak_Promoter, 47
K562, GRO-cap, 3_Poised_Promoter, 7
K562, dREG, 3_Poised_Promoter, 12
K562, GRO-cap, 4_Strong_Enhancer, 52
K562, dREG, 4_Strong_Enhancer, 65
K562, GRO-cap, 5_Strong_Enhancer, 13
K562, dREG, 5_Strong_Enhancer, 28
K562, GRO-cap, 6_Weak_Enhancer, 9
K562, dREG, 6_Weak_Enhancer, 20
K562, GRO-cap, 7_Weak_Enhancer, 2
K562, dREG, 7_Weak_Enhancer, 7
GM12878, GRO-cap, 1_Active_Promoter, 82
GM12878, dREG, 1_Active_Promoter, 87
GM12878, GRO-cap, 2_Weak_Promoter, 17
GM12878, dREG, 2_Weak_Promoter, 39
GM12878, GRO-cap, 3_Poised_Promoter, 9
GM12878, dREG, 3_Poised_Promoter, 16
GM12878, GRO-cap, 4_Strong_Enhancer, 46
GM12878, dREG, 4_Strong_Enhancer, 61
GM12878, GRO-cap, 5_Strong_Enhancer, 11
GM12878, dREG, 5_Strong_Enhancer, 27
GM12878, GRO-cap, 6_Weak_Enhancer, 6
GM12878, dREG, 6_Weak_Enhancer, 15
GM12878, GRO-cap, 7_Weak_Enhancer, 1
GM12878, dREG, 7_Weak_Enhancer, 4", header=TRUE)

library(RColorBrewer)

myColours <- brewer.pal(6,"Blues")

my.settings <- list(
  superpose.polygon=list(col=myColours[c(5,2)], border="transparent"),
  strip.background=list(col=myColours[6]),
  strip.border=list(col="black")
)

barchart(Sensitivity~Ernst.classes | Cell.type, groups=Method, data=acc, ylim=c(0,100), 
		scales=list(x=list(rot=65)), par.settings = my.settings, par.strip.text=list(col="white", font=2),
		auto.key=list(space="top", columns=2, title="Experiment", cex.title=1))

dev.off()
