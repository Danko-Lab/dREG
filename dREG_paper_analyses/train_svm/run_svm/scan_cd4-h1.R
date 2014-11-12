require(featureDetector)

## Read PRO-seq data.
gs_plus_path  <- c("J-U-cc_plus.bw") #c("H1-U_plus.bw", "J-U_plus.bw", "H1-PI_plus.bw", "J-PI_plus.bw")
gs_minus_path <- c("J-U-cc_minus.bw") #c("H1-U_minus.bw", "J-U_minus.bw", "H1-PI_minus.bw", "J-PI_minus.bw")

outnames <- c("J-U.predictions.bedGraph") #c("H1-U.predictions.bedGraph", "J-U.predictions.bedGraph", "H1-PI.predictions.bedGraph", "J-PI.predictions.bedGraph")

load("asvm.RData")
gdm <- genomic_data_model(window_sizes= c(10, 25, 50, 500, 5000), half_nWindows= c(10, 10, 30, 20, 20)) 

for(i in 1:length(gs_plus_path)) {
	## Now scan all positions in the genome ...
	inf_positions <- get_informative_positions(gs_plus_path[i], gs_minus_path[i], depth= 0, step=50, use_ANDOR=TRUE, use_OR=FALSE) ## Get informative positions.

	pred_val<- eval_reg_svm(gdm, asvm, inf_positions, gs_plus_path[i], gs_minus_path[i], batch_size= 10000, ncores=32)

	final_data <- data.frame(inf_positions, pred_val)
	options("scipen"=100, "digits"=4)
	write.table(final_data, file=outnames[i], row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
}
