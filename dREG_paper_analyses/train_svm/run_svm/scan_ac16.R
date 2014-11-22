require(dREG)

load("asvm.intersDNase.getTrainSet.RData")

## Read PRO-seq data.
ps_plus_path  <- "ac16.unt.all_plus.bw"
ps_minus_path <- "ac16.unt.all_minus.bw"

## Now scan all positions in the genome ...
inf_positions <- get_informative_positions(ps_plus_path, ps_minus_path, depth= 0, step=50, use_ANDOR=TRUE, use_OR=FALSE) ## Get informative positions.
gdm <- genomic_data_model(window_sizes= c(10, 25, 50, 500, 5000), half_nWindows= c(10, 10, 30, 20, 20)) 

pred_val<- eval_reg_svm(gdm, asvm, inf_positions, ps_plus_path, ps_minus_path, batch_size= 50000, ncores=62)

final_data <- data.frame(inf_positions, pred_val)
options("scipen"=100, "digits"=4)
write.table(final_data, file="ac16.predictions.bedGraph", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")


