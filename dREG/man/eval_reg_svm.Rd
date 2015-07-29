\name{eval_reg_svm}
\alias{eval_reg_svm}
\title{ Evaluates a set of genomic coordinates for regulatory potential using P/GRO-seq data }
\description{ Evaluates a set of genomic coordinates for regulatory potential using P/GRO-seq data }
\usage{
eval_reg_svm(gdm, 
      asvm, 
      positions, 
      bw_plus_path, 
      bw_minus_path, 
      batch_size = 50000, 
      ncores = 3, 
      debug = TRUE)
}
\arguments{
  \item{gdm}{Genomic data model return by \code{\link{genomic_data_model}}. }
  \item{asvm}{A pre-trained SVM model from the e1071 package returned by \code{\link{regulatory_svm}}.}
  \item{positions}{Data frame with 2 columns indicating the universe of positions to test and evaluate(chrom,chromCenter). It can be returned by \code{\link{get_informative_positions}}.}
  \item{bw_plus_path}{String value indicating file path to bigWig file representing the plus strand.}
  \item{bw_minus_path}{String value indicating file path to bigWig file representing the minus strand}
  \item{batch_size}{Number of positions to evaluate at once (more might be faster, but takes more memory).}
  \item{ncores}{The number of cpu cores in parallel computing}
  \item{debug}{Logical value indicating the process detail is outputted.}
}
\details{
}
\value{
Returns the value of the SVM for each genomic coordinate specified.
}
\references{
Danko, C. G., Hyland, S. L., Core, L. J., Martins, A. L., Waters, C. T., Lee, H. W., ... & Siepel, A. (2015). 
Identification of active transcriptional regulatory elements from GRO-seq data. Nature methods, 12(5), 433-438.
}
\seealso{
\code{\link{get_informative_positions}}, \code{\link{get_test_set}}, \code{\link{read_genomic_data}}, \code{\link{regulatory_svm}}
}
\examples{


## The following codes cannot run without the bigWig files

# ps_plus_path  <- "bigwig.plus.bw"
# ps_minus_path <- "bigwig.minus.bw"

## Now scan all positions in the genome ...
# positions <- get_informative_positions(ps_plus_path, ps_minus_path, 
#		depth= 0, step=50, use_ANDOR=TRUE, use_OR=FALSE);

# pred_val<- eval_reg_svm( gdm, asvm, inf_positions, ps_plus_path, ps_minus_path, batch_size= 50000)
# write.table( data.frame(inf_positions, pred_val), file="eval.tab", 
#	row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")


}
\keyword{ function }
