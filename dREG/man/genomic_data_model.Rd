\name{genomic_data_model}
\alias{genomic_data_model}
\title{
Creates a genome data model.
}
\description{
Creates a genome data model.
}
\usage{
genomic_data_model(window_sizes, half_nWindows)
}
\arguments{
  \item{window_sizes}{}
  \item{half_nWindows}{}
}
\details{
}
\value{
A s4 object is returned with 
\item{n_zooms}{Number indicating zoom ratio.}
\item{window_sizes}{Vector indicating window sizes.}
\item{half_nWindows}{Vector indicating number of half windows.}
}
\references{
Danko, C. G., Hyland, S. L., Core, L. J., Martins, A. L., Waters, C. T., Lee, H. W., ... & Siepel, A. (2015). 
Identification of active transcriptional regulatory elements from GRO-seq data. Nature methods, 12(5), 433-438.
}
\seealso{
\code{\link{read_genomic_data}}, \code{\link{regulatory_svm}},  \code{\link{eval_reg_svm}}
}
\examples{
gdm <- genomic_data_model( c(10,20,30), c(10, 10, 10) )
}
\keyword{ function }

