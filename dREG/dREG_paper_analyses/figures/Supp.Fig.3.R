
pdf("figS3.pieCharts.pdf")

	slices <- c(1936, 6503, 11336, 8952, 1339, 640, 19287, 4343)
	lbls <- c("10_Txn_Elongation", "11_Weak_Txn", "12_Repressed", "13_Heterochrom/lo", "14_Repetitive/CNV", "15_Repetitive/CNV", "8_Insulator", "9_Txn_Transition")
	pie(slices, labels = lbls, main="Functional elements associated with DNAse-1-unique regions")

	slices <- c(sum(1936, 6503, 4343), sum(11336, 8952), sum(1339, 640), 19287)
	lbls <- c("Transcribed", "Repressed/ Heterochrom", "Repetitive", "Insulator")
	pie(slices, labels = lbls, main="Functional elements associated with DNAse-1-unique regions")


	slices <- c(904, 4321, 2785, 2129, 5126, 17814, 56986)
	lbls <- c("1_Active_Promoter", "2_Weak_Promoter", "3_Poised_Promoter", "4_Strong_Enhancer", "5_Strong_Enhancer", "6_Weak_Enhancer", "7_Weak_Enhancer")
	pie(slices, labels = lbls, main="chromHMM regions without DNAse-1 hypersensitivity")

	slices <- c(sum(904, 4321, 2785), sum(2129, 5126), sum(17814, 56986))
	lbls <- c("Promoter", "Strong_Enhancer", "Weak_Enhancer")
	pie(slices, labels = lbls, main="chromHMM regions without DNAse-1 hypersensitivity")

	slices <- c(sum(904, 4321, 2785), sum(2129, 5126), sum(17814, 56986))
	lbls <- c("Promoter", "Strong_Enhancer", "Weak_Enhancer")
	pie(slices, labels = lbls, main="chromHMM regions without DNAse-1 hypersensitivity")


	slices <- c(sum(904, 4321, 2785), sum(2129, 5126), sum(17814, 56986))
	lbls <- c("Promoter", "Strong_Enhancer", "Weak_Enhancer")
	pie(slices, labels = lbls, main="chromHMM regions without DNAse-1 hypersensitivity")


	slices <- c(4221, 64806)
	lbls <- c("Other", "Weak_Enhancer")
	pie(slices, labels = lbls, main="chromHMM regions without DNAse-1 hypersensitivity")

dev.off()

