\name{roc.calc}
\alias{roc.calc}
\title{
Calculates the TPR and FPR for a ROC plot.
}
\description{
Calculates the TPR and FPR for a ROC plot.
}
\usage{
roc.calc(true, 
      possible, 
      scores, 
      filterPossible = TRUE, 
      n_points = 100)
}
\arguments{
  \item{true}{Bed-style data frame, a set of 'true' genomic intervals (e.g. ChIP-seq peaks).}
  \item{possible}{Bed-style data frame, A set of 'possible' genomic intervals (e.g. DNAse-1 peaks).}
  \item{scores}{Vector indicating the scores for each possibe genomic interval in parameter of \code{possible}. }
  \item{filterPossible}{Vector indicating indexes which be removed.}
  \item{n_points}{Integer indicating how many points for the ROC plot.}
}
\details{
}
\value{
A data frame with 3 columns is returned
\item{FPR}{False Positive Rate}
\item{TPR}{True Positive Rate}
\item{threshold}{Threshold based on the score parameter.}

}
\references{
Danko, C. G., Hyland, S. L., Core, L. J., Martins, A. L., Waters, C. T., Lee, H. W., ... & Siepel, A. (2015). 
Identification of active transcriptional regulatory elements from GRO-seq data. Nature methods, 12(5), 433-438.
}
\seealso{
\code{\link{logreg.roc.calc}}, \code{\link{combine.roc}}, \code{\link{roc.auc}}, \code{\link{roc.plot}}
}
\examples{
}
\keyword{ plot }
