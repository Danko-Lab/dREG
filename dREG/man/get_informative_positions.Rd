\name{get_informative_positions}
\alias{get_informative_positions}
\title{Gets center positions that pass a minimum depth filter}
\description{Returns a data frame with center positions that pass a minimum depth filter}
\usage{
get_informative_positions(bw_path, 
      bw_minus_path = NULL, 
      depth = 0, 
      window = 400, 
      step = 50, 
      use_OR = TRUE, 
      use_ANDOR = TRUE, 
      debug = TRUE)
}
\arguments{
  \item{bw_path}{String indicating file path to bigwig file representing the plus strand.}
  \item{bw_minus_path}{String indicating file path to bigwig file representing the minus strand, If specified, takes the windows that pass the step in both bigWig files.(intersection)}
  \item{depth}{Integer value indicating minimum number of reads to return.}
  \item{window}{Integer value indicating window distance between to search for #depth reads [bp].}
  \item{step}{Integer value indicating step distance for window list.}
  \item{use_OR}{Logical value indicating if the center positions in minus bigwig file are merged into the results. If false, the intersection operation will be performed to the center positions of plus bigwig and from minus bigwig.}
  \item{use_ANDOR}{Logical value indicating if the center positions will be merged from the two results.  a) Intersection operation with the conditions: window interval=1000 depth>=0. b) Union operation with with the conditions: window interval=100 depth >=2.}
  \item{debug}{Logical value indication the process detail is outputted.}
}
\details{
The use_ANDOR and use_OR parameter are applied to two Bigwig files as following logical:\cr

\preformatted{
if(use_ANDOR){
    v1 <- get_window_Or  (window=1000, depth=0);
    v2 <- get_window_and (window=100,  depth=2);
    vals <- c(v1,v2);
}
else {
  if(use_OR){
     vals <- get_window_Or( window=window, depth=depth);
  }
  else {
     vals <- get_window_and( window=window, depth=depth);
  }
}
}

}

\value{
A BED-style data frame will be returned with 3 columns
\item{chrom }{Chromosome information}
\item{chromStart }{Start position}
\item{chromEnds }{End position}

}
\references{
Danko, C. G., Hyland, S. L., Core, L. J., Martins, A. L., Waters, C. T., Lee, H. W., ... & Siepel, A. (2015). 
Identification of active transcriptional regulatory elements from GRO-seq data. Nature methods, 12(5), 433-438.
}

\seealso{
\code{\link{get_test_set}}, \code{\link{read_genomic_data}}, \code{\link{regulatory_svm}}, \code{\link{eval_reg_svm}}
}
\keyword{ function }