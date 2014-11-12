##
## deep_generative_model.R -- 
## Trains a generative DBN to detect distinct regulatory region types.
##
##

nsamp <- 20000

shuffle10 <- function(n_elements) {
  indx <- c(1:n_elements)
  shuf <- c(which(indx %% 10 == 9), which(indx %% 10 == 8), which(indx %% 10 == 7), which(indx %% 10 == 6), which(indx %% 10 == 5), 
			which(indx %% 10 == 4), which(indx %% 10 == 3), which(indx %% 10 == 2), which(indx %% 10 == 1), which(indx %% 10 == 0))
  return(order(shuf))
}

#setwd("/usr/projects/GROseq.parser/tss_detecter/")
setwd("~/work/dbn_shapes")
require(featureDetector)

## Read PRO-seq data & positive regions.
## Read PRO-seq data.
 ps_plus_path  <- "/usr/data/GROseq.parser/hg19/k562/proseq/K562_unt.sort.bed.gz_plus.bw" #K562_unt.subsamp10pct.bed.gz_plus.bw"
 ps_minus_path <- "/usr/data/GROseq.parser/hg19/k562/proseq/K562_unt.sort.bed.gz_minus.bw" #K562_unt.subsamp10pct.bed.gz_minus.bw"

## Get informative positions for training.
 inf_positions <- get_informative_positions(ps_plus_path, ps_minus_path, depth= 0, window= 500, step=100, use_OR=TRUE) ## Get informative positions.
 print(paste("Number of inf. positions: ", NROW(inf_positions)))

 ## Intersect with GRO-cap to get start sites.
 tss <- read.table("training_beds/grocaptss.bed")

 ## And with chromHMM to get active gene bodies.
 txn <- read.table("training_beds/txn.strand.bed")

 ## And with chromHMM to get active gene bodies.
 ends <- read.table("training_beds/geneEnd.bed")

 
## Search the space of possible representations.
 gdm <- genomic_data_model(window_sizes= c(10, 25, 50, 500, 5000), half_nWindows= c(10, 10, 30, 20, 20))
 adbn <- dbn(layer_sizes= c(360,1000,500,50), batch_size=10, cd_n=1, momentum_decay= 0.9, weight_cost= 1e-5, learning_rate=0.1)
 indx <- shuffle10(2*nsamp)

## Cut down here a little bit...
set.seed(61)
for(i in 1:100) {
 ip_sample <- rbind(get_test_set(inf_positions, tss, nsamp), get_test_set(inf_positions, txn[txn$V6=="+",], nsamp),  get_test_set(inf_positions, txn[txn$V6=="-",], nsamp),
					get_test_set(inf_positions, ends[ends$V6=="+",], nsamp),  get_test_set(inf_positions, ends[ends$V6=="-",], nsamp))[indx,]

 data_mat <- read_genomic_data(gdm, ip_sample, ps_plus_path, ps_minus_path)
 adbn <- dbn.pretrain(adbn, data= t(data_mat), n_epocs= 1, n_threads=8)
}

## Refine.
indx_to_refine <- 1:(2*nsamp)
refine_sample <- rbind(get_test_set(inf_positions, tss, 2*nsamp)[indx_to_refine,], 
			get_test_set(inf_positions, txn[txn$V6=="+",], 2*nsamp)[indx_to_refine,],  get_test_set(inf_positions, txn[txn$V6=="-",], 2*nsamp)[indx_to_refine,],
			get_test_set(inf_positions, ends[ends$V6=="+",], 2*nsamp)[indx_to_refine,],  get_test_set(inf_positions, ends[ends$V6=="-",], 2*nsamp)[indx_to_refine,])[indx,]
cor_labels <- as.factor(c(rep("tss", 2*nsamp), rep("txn_plus", 2*nsamp), rep("txn_minus", 2*nsamp), rep("end_plus", 2*nsamp), rep("end_minus", 2*nsamp))[indx])
					
rdbn <- dbn.refine(adbn, data= refine_sample, labels=cor_labels, n_epocs=10)

save.image("~/work/dbn_shapes/proseq.adbn.forDBN.1.RData")

## Then run on a small test set!!!



#
#pdf("proseq.model.receptive_fields.pdf")
#for(i in 1:18) {
#  plot.gdm(gdm, gdata= t(dbn.clamplayer(adbn, neuron= i, layer= 5)))
#}
#dev.off()

