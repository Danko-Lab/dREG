## train_svm.R -- Trains an SVM to detect regulatory regions, and uses it to.
## scan the genome in 10bp windows.
##
args<-commandArgs(TRUE)
require(featureDetector)
filename <- args[3]             # this will become the filename for a pdf later

## Read PRO-seq data.
ps_plus_path <- args[1]         # these should be bws
ps_minus_path <- args[2]

## Get positive regions.
GROcap_tss_bed <- read.table("/usr/projects/GROseq.parser/tss_detecter/andre_hmm/hg19.k562.new_hmm2b.post2.bed") #"/usr/projects/GROseq.parser/tss_new/hg19.k562.new_hmm2.bed", skip=1)

## Train the SVM.
inf_positions <- get_informative_positions(ps_plus_path, ps_minus_path, depth= 0, step=50, use_ANDOR=TRUE, use_OR=FALSE) ## Get informative positions.
print(paste("Number of inf. positions: ", NROW(inf_positions)))

gdm <- genomic_data_model(window_sizes= c(10, 25, 50, 500, 5000), half_nWindows= c(10, 10, 30, 20, 20)) 
extra_enrich_bed <- read.table("/usr/projects/GROseq.parser/tss_detecter/GencodeMerge.IntersectOpStrand.bed")
allow_bed <- read.table("/usr/projects/GROseq.parser/tss_detecter/chromHmm.k562.enh.prom.bed.gz")
asvm <- regulatory_svm(gdm, ps_plus_path, ps_minus_path, inf_positions, GROcap_tss_bed, pdf_path= "roc_plot.and1.lgModel.pdf", n_train=25000, n_eval=1000, extra_enrich_bed= extra_enrich_bed, allow= allow_bed)

## Now scan all positions in the genome ...
inf_positions <- get_informative_positions(ps_plus_path, ps_minus_path, depth= 0, step=50, use_ANDOR=TRUE, use_OR=FALSE) ## Get informative positions.

pred_val<- eval_reg_svm(gdm, asvm, inf_positions, ps_plus_path, ps_minus_path, batch_size= 10000, ncores=16)

final_data <- data.frame(inf_positions, pred_val)
options("scipen"=100, "digits"=4)
write.table(final_data, file="predictions.bedGraph", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

save.image(paste(filename,".asvm.RData",sep=""))

###################################
## Analyze the trained model ... 

##  Finally, we'll merge consecutive positive values and intersect with Andre's bigWig.
#final_data <- read.table("predictions.bedGraph")           # don't need to do this, already have this information
thre <- 0.7
#final_data <- final_data[final_data[,4]> thre,]
pred_pos <- final_data[final_data[,4]> thre,]
write.table(pred_pos, file="predictions.bed", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

system("cat predictions.bed | awk 'BEGIN{OFS=\"\t\"} {print $1,$2-151,$3+151,\"N\",$4'} | bedops --merge - > predictions.merge.bed", intern=TRUE)
system("cat predictions.bedGraph | awk 'BEGIN{OFS=\"\t\"} {print $1,$2-50,$3+50,\"N\",$4'} | bedmap  --echo --max predictions.merge.bed - | sed s/\\|/\"\t\"/g > predictions.bed", intern=TRUE)

final_pred_pos <- read.table("predictions.bed")

GROcap_tss_bed <- read.table("/usr/projects/GROseq.parser/tss_new/hg19.k562.new_hmm2.bed",skip=1) #read.table("/usr/projects/GROseq.parser/tss_new/hg19.k562.new_hmm2.bed", skip=1)

rocPlot <- roc.calc(GROcap_tss_bed, final_pred_pos, final_pred_pos[,4], filterPossible=TRUE)

desired_fpr<-0.1
the_index<-which(abs(rocPlot[,1]-desired_fpr)==min(abs(rocPlot[,1]-desired_fpr)))
tpr<-rocPlot[the_index,2]
fpr<-rocPlot[the_index,1]
auc<-roc.auc(rocPlot)

results<-data.frame(filename,tpr,fpr,auc)
write.table(results,file='sens_results.txt',row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE,append=TRUE)
print(paste(fpr,tpr,roc.auc(rocPlot)))

pdf(paste(filename,'.pdf',sep=''))
roc.plot(rocPlot)
dev.off()
