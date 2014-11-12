require(bigWig)
require(featureDetector)

bw_plus  <- "M3-U.bed.gz_plus.bw" 
bw_minus <- "M3-U.bed.gz_minus.bw"

inf_positions <- get_informative_positions(bw_plus, bw_minus)

