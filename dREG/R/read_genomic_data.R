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
read_genomic_data <- function(gdm, bed, bigwig_plus, bigwig_minus, as_matrix= TRUE) {

  stopifnot(NROW(gdm@window_sizes) == NROW(gdm@half_nWindows))
  zoom<- list(as.integer(gdm@window_sizes), as.integer(gdm@half_nWindows))

  dat <- .Call("get_genomic_data_R", as.character(bed[,1]), as.integer(floor((bed[,3]+bed[,2])/2)), as.character(bigwig_plus), as.character(bigwig_minus), zoom, PACKAGE= "dREG")

  if(as_matrix) 
    dat <- t(matrix(unlist(dat), ncol=NROW(bed)))
    #dat <- t(matrix(unlist(lapply(c(1:NROW(dat)), function(x) {unlist(dat[[x]])})), ncol=NROW(bed)))

	return(dat)
}


