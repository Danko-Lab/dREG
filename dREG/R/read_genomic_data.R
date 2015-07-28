## read_genomic_data -- reads genomic data in the specified zoom parameters...
##

setClass("genomic_data_model",#"restricted_boltzman_machine", 
  representation(
    n_zooms="integer",
    window_sizes="integer",
    half_nWindows="integer"
    ),
)

genomic_data_model <- function(window_sizes, half_nWindows) {
  stopifnot(NROW(window_sizes) == NROW(half_nWindows))
  new("genomic_data_model", n_zooms= as.integer(NROW(window_sizes)), window_sizes= as.integer(window_sizes), half_nWindows= as.integer(half_nWindows))
}

#' Reads genomic data from the specified position...
#'
#' @param bed A data.frame of genomic regions in bed format (chrom, start, end).
#' @param bigwig_plus Path to bigwig file representing GRO-seq/ PRO-seq reads on the plus strand.
#' @param bigwig_minus Path to bigwig file representing GRO-seq/ PRO-seq reads on the minus strand.
#' @param as_matrix If true, returns a matrix object.
#' @return Returns a list() object, where each element in the list is the zoom data
#' centered on a 
read_genomic_data <- function(gdm, bed, file_bigwig_plus, file_bigwig_minus, as_matrix= TRUE, scale.method=c("logistic", "totalRead")) {
    if(missing(scale.method)){scale.method<-"logistic"};
  stopifnot(NROW(gdm@window_sizes) == NROW(gdm@half_nWindows))
  zoom<- list(as.integer(gdm@window_sizes), as.integer(gdm@half_nWindows))
  if(scale.method=="logistic"){
  dat <- .Call("get_genomic_data_R", as.character(bed[,1]), as.integer(floor((bed[,3]+bed[,2])/2)), as.character(file_bigwig_plus), as.character(file_bigwig_minus), zoom, as.logical(TRUE), PACKAGE= "dREG")
  }
  else{
      dat <- .Call("get_genomic_data_R", as.character(bed[,1]), as.integer(floor((bed[,3]+bed[,2])/2)), as.character(file_bigwig_plus), as.character(file_bigwig_minus), zoom, as.logical(FALSE), PACKAGE= "dREG")
      total.read.count<- sum(abs(get_reads_from_bigwig(file_bigwig_plus, file_bigwig_minus)));
      dat<-lapply(dat, "/", total.read.count);
  }
  
  if(as_matrix) 
    dat <- t(matrix(unlist(dat), ncol=NROW(bed)))
    
    #dat <- t(matrix(unlist(lapply(c(1:NROW(dat)), function(x) {unlist(dat[[x]])})), ncol=NROW(bed)))

	return(dat)
}

# query read counts of all chromosomes from bigWig file.
#
# @bw.plus,     bigWig object
# @bw.minus,    bigWig object
# @chromInfo    data.frame with 2 columns(chr, size);
#
# @return       vector of reads in plus and minus file.
get_reads_from_bigwig <- function(file_bigwig_plus, file_bigwig_minus)
{
    bw.plus<-load.bigWig(file_bigwig_plus)
    bw.minus<-load.bigWig(file_bigwig_minus)
    offset_dist <- 250;
    df.bed.plus<-data.frame(bw.plus$chroms, offset_dist, bw.plus$chromSizes, names=".", scores=".",strands="+")
    df.bed.minus<-data.frame(bw.minus$chroms, offset_dist, bw.minus$chromSizes, names=".", scores=".", strands="-")
    r.plus <- sum(abs(bed6.region.bpQuery.bigWig( bw.plus, bw.minus, df.bed.plus)));
    r.minus <- sum(abs(bed6.region.bpQuery.bigWig( bw.plus, bw.minus, df.bed.minus)));
        
    return(c(r.plus,r.minus));
}