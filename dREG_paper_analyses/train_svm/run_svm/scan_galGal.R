require(featureDetector)

## Read PRO-seq data.
gs_plus_path  <- c("2106_5035_6611_N_chEmb-1_WE_TTAGGC_R1.pl.bw", "2106_5035_6612_N_chEmb-2_HD_TGACCA_R1.pl.bw", "2106_5035_6613_N_chEmb-3_LH_ACAGTG_R1.pl.bw", "2106_5035_6614_N_chEmb-4_RH_GCCAAT_R1.pl.bw", "2106_5035_6615_N_chEmb-5_LM_CAGATC_R1.pl.bw", "2106_5035_6616_N_chEmb-6_RM_ACTTGA_R1.pl.bw")
gs_minus_path <- c("2106_5035_6611_N_chEmb-1_WE_TTAGGC_R1.mn.bw", "2106_5035_6612_N_chEmb-2_HD_TGACCA_R1.mn.bw", "2106_5035_6613_N_chEmb-3_LH_ACAGTG_R1.mn.bw", "2106_5035_6614_N_chEmb-4_RH_GCCAAT_R1.mn.bw", "2106_5035_6615_N_chEmb-5_LM_CAGATC_R1.mn.bw", "2106_5035_6616_N_chEmb-6_RM_ACTTGA_R1.mn.bw")

outnames <- c("chEmb-1_WE_TTAGGC.TSS.bedGraph", "chEmb-2_HD_TGACCA.TSS.bedGraph", "chEmb-3_LH_ACAGTG.TSS.bedGraph", "chEmb-4_RH_GCCAAT.TSS.bedGraph", "chEmb-5_LM_CAGATC.TSS.bedGraph", "chEmb-6_RM_ACTTGA.TSS.bedGraph")

load("asvm.RData")
gdm <- genomic_data_model(window_sizes= c(10, 25, 50, 500, 5000), half_nWindows= c(10, 10, 30, 20, 20)) 

for(i in 1:length(gs_plus_path)) {
	## Now scan all positions in the genome ...
	inf_positions <- get_informative_positions(gs_plus_path[i], gs_minus_path[i], depth= 0, step=50, use_ANDOR=TRUE, use_OR=FALSE) ## Get informative positions.

	pred_val<- eval_reg_svm(gdm, asvm, inf_positions, gs_plus_path[i], gs_minus_path[i], batch_size= 10000, ncores=62)

	final_data <- data.frame(inf_positions, pred_val)
	options("scipen"=100, "digits"=4)
	write.table(final_data, file=outnames[i], row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
}
