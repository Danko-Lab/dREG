\name{read_genomic_data}
\alias{read_genomic_data}
\title{Gets read data from the specified genomic position.}
\description{Gets read data from the specified genomic position.}
\usage{
read_genomic_data( gdm, 
      bed, 
      file_bigwig_plus, 
      file_bigwig_minus, 
      as_matrix = TRUE, 
      scale.method = c("logistic", "linear"),
      batch_size = 50000, 
      ncores = 1 )
}
\arguments{
  \item{gdm}{Genomic data model return by \code{\link{genomic_data_model}}.}
  \item{bed}{bed-style data frame of genomic regions.(at least 3 columns including chrom, start, end).}
  \item{file_bigwig_plus}{String value indicating file path to bigwig file representing GRO-seq/ PRO-seq reads on the plus strand.}
  \item{file_bigwig_minus}{String value indicating file path to bigwig file representing GRO-seq/ PRO-seq reads on the minus strand.}
  \item{as_matrix}{Logical type,if true, returns a matrix object, otherwise returns a list() object, where each element in the list is the zoom data.}
  \item{scale.method}{String value indicating the normalize method of read counts. Two options are available, "logistic" or "linear", default value is logistic. See details}
  \item{batch_size}{Number of genomic positions to evaluate at once (more might be faster, but takes more memory)}
  \item{ncores}{Number of CPU cores in parallel computing}
}
\details{
Data normalize method: \cr
(1): Logistic function: F(x) = 1/(1+exp(-a*(x-b))\cr
(2): Linear function: F(x) = x / tootal_reads\cr
}
\value{
A matrix of normalized read count, the columns are windows list specified by \code{gdm} object.
}
\references{
Danko, C. G., Hyland, S. L., Core, L. J., Martins, A. L., Waters, C. T., Lee, H. W., ... & Siepel, A. (2015). 
Identification of active transcriptional regulatory elements from GRO-seq data. Nature methods, 12(5), 433-438.
}
\seealso{
\code{\link{get_informative_positions}}, \code{\link{get_test_set}}, \code{\link{regulatory_svm}}, \code{\link{eval_reg_svm}}
}
\examples{

file_bigwig_plus <- "";
file_bigwig_minus <- "";
gdm <- genomic_data_model(20, 10);
#mat <- read_genomic_data(gdm, bed, file_bigwig_plus, file_bigwig_minus);
#summary(mat);

}
\keyword{ function }
