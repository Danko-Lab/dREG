\name{roc.plot}
\alias{roc.plot}
\title{Draws a ROC figure.
}
\description{Draws a ROC figure.}
\usage{
roc.plot(ROC, ...)
}
\arguments{
  \item{ROC}{Matrix or data frame with 3 columns, FPR, TPR and threshold.}
  \item{\dots}{The parameters for \code{plot} function}
}
\value{
None
}
\references{
Danko, C. G., Hyland, S. L., Core, L. J., Martins, A. L., Waters, C. T., Lee, H. W., ... & Siepel, A. (2015). 
Identification of active transcriptional regulatory elements from GRO-seq data. Nature methods, 12(5), 433-438.
}

\seealso{
\code{\link{roc.calc}}, \code{\link{logreg.roc.calc}}, \code{\link{combine.roc}}, \code{\link{roc.auc}}
}
\examples{

true <- c(rep(1, 100), rep(0, 100));
scores <- c( rnorm(100, 1, 1 ), rnorm(100, 0, 1 ) );
roc_mat <- logreg.roc.calc( true, scores );
AUC<- roc.auc(roc_mat);
roc.plot(roc_mat, main=AUC );

}
\keyword{ plot }
