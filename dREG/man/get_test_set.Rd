\name{get_test_set}
\alias{get_test_set}
\title{ Returns a genome loci of positive set and negative set for SVM training purpose.}
\description{ Returns a genome loci of positive set and negative set for SVM training purpose }
\usage{
get_test_set(positions, 
      positive, 
      n_samp, 
      allow = NULL, 
      enrich_negative_near_pos = 0.15, 
      extra_enrich_bed = NULL, 
      extra_enrich_frac = 0.1, 
      avoid_dist = 100)
}
\arguments{
  \item{positions}{Bed-style data frame indicating the universe of positions to test and evaluate (chrom,chromCenter).}
  \item{positive}{Bed-style data frame containing positive positions (chrom,chromStart,chromEnd).}
  \item{n_samp}{Number of training examples}
  \item{allow}{Bed-style data frame containing inverse negative set of positions (chrom,chromStart,chromEnd).}
  \item{enrich_negative_near_pos}{Fraction of training examples chosen to be nearby (<=5kb) a positive example [0,1].}
  \item{extra_enrich_bed}{Bed-style data frame indicating extra bed file to enrich near.}
  \item{extra_enrich_frac}{Fraction of final positions sampled in the negative set which are in the bed file. Unused if extra_enrich_bed is NULL.}
  \item{avoid_dist}{Integer value indicating how long extend avoiding genomic loci.}
}
\details{
(1). The parameter of \code{positions} can be obtained by \code{\link{get_informative_positions}}.

}
\value{
Returns a data frame including double number of the _train set(2*n_samp), each sample includes 4 items.
  \item{chrom}{}
  \item{chromStart}{}
  \item{chromEnd}{}
  \item{status}{1 for positive and 0 for negative.}

}
\references{
Danko, C. G., Hyland, S. L., Core, L. J., Martins, A. L., Waters, C. T., Lee, H. W., ... & Siepel, A. (2015). 
Identification of active transcriptional regulatory elements from GRO-seq data. Nature methods, 12(5), 433-438.
}
\seealso{
\code{\link{get_informative_positions}}, \code{\link{read_genomic_data}}, \code{\link{regulatory_svm}}, \code{\link{eval_reg_svm}}

}

\keyword{ function }
