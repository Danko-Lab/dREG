##
## Test combo R/C function
#GROcap_tss_bed <- read.table("/usr/projects/GROseq.parser/tss_new/hg19.k562.new_hmm2.bed", skip=1)
require(Rdbn)

gs_plus  <- "/usr/data/GROseq.parser/hg19/k562/proseq_celastrol_prelim/celastrol_proseq_0min_plus.bigWig"
gs_minus <- "/usr/data/GROseq.parser/hg19/k562/proseq_celastrol_prelim/celastrol_proseq_0min_minus.bigWig"
bw_data_plus <- load.bigWig(gs_plus)
bw_data_minus <- load.bigWig(gs_minus)

step_size=10
n_windows=20

Gencode <- read.table("/usr/projects/GROseq.parser/annotations/gencode.comprehensive.bed", header=FALSE, skip=1)
Gencode <- Gencode[Gencode[,11]== "protein_coding",c(1:3)]
Gencode <- Gencode[sample(c(1:NROW(Gencode)), 10000),]

bigWig_gencode_data <- colSums(collect.many(Gencode, bw_data_plus, bw_data_minus, halfWindow= step_size*n_windows, step= step_size))
Rdbn_gencode_data <- colSums(read_genomic_data(Gencode, gs_plus, gs_minus, window_sizes= step_size, half_nWindows= n_windows))

plot(bigWig_gencode_data, type="l", ylim=c(0, max(c(bigWig_gencode_data, Rdbn_gencode_data))))
points(Rdbn_gencode_data, type="l", col="red")



## Found a bug when one of the windows overlaps the 0bp chromosome boundary.  This code (below) is sufficient to reproduce it.
## UPDATE: bug traced to bigwig_read_i not returning data when start<0.  Should now be fixed.

 require(featureDetector)
 ps_plus_path  <- "/usr/data/GROseq.parser/hg19/k562/proseq/K562_unt.sort.bed.gz_plus.bw" #/usr/data/GROseq.parser/hg19/k562/proseq/
 ps_minus_path <- "/usr/data/GROseq.parser/hg19/k562/proseq/K562_unt.sort.bed.gz_minus.bw" #/usr/data/GROseq.parser/hg19/k562/proseq/
 
 gdm <- genomic_data_model(window_sizes= c(10, 10000), half_nWindows= c(10, 2)) 
 gdm2 <- genomic_data_model(window_sizes= c(10, 10000), half_nWindows= c(10, 1)) 

 inf_positions <- read.table("tmp.bedgraph")
 print(paste("Number of inf. positions: ", NROW(inf_positions)))
 
 colSums(read_genomic_data(gdm, inf_positions, ps_plus_path, ps_minus_path, as_matrix= TRUE))
 colSums(read_genomic_data(gdm2, inf_positions, ps_plus_path, ps_minus_path, as_matrix= TRUE))
 
 plot(colSums(read_genomic_data(gdm, inf_positions, ps_plus_path, ps_minus_path, as_matrix= TRUE)), type="l")
 points(colSums(read_genomic_data(gdm2, inf_positions, ps_plus_path, ps_minus_path, as_matrix= TRUE)), type="l", col="gray")

 ## Actual per base values from the bigWig.
 bwp <- load.bigWig(ps_plus_path)
 bwm <- load.bigWig(ps_minus_path)
 query.bigWig(bwp, "chr1", 9100, 11100)
 query.bigWig(bwm, "chr1", 9100, 11100)
 
 