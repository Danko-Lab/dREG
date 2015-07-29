\name{logreg.roc.calc}
\alias{logreg.roc.calc}
\title{
Calculates the TPR and FPR for a ROC plot.
}
\description{
Calculates the TPR and FPR for a ROC plot from the status and score vector.
}
\usage{
logreg.roc.calc(true, scores)
}
\arguments{
  \item{true}{Vector indicating the two status, 1 and 0. }
  \item{scores}{Vector indicating the scores for each status calculated by the predict function.}
}
\details{

The function of \code{\link{roc.calc}} calculates a ROC matrix for the genomic loci, 
whereas the function of \code{\link{logreg.roc.calc}} calculates for a status vector.
}
\value{
A data frame with 3 columns is returned, which is same as \code{\link{roc.calc}}.
\item{FPR}{False Positive Rate.}
\item{TPR}{True Positive Rate.}
\item{threshold}{Threshold based on the score parameter.}
}
\references{
Danko, C. G., Hyland, S. L., Core, L. J., Martins, A. L., Waters, C. T., Lee, H. W., ... & Siepel, A. (2015). 
Identification of active transcriptional regulatory elements from GRO-seq data. Nature methods, 12(5), 433-438.
}
\seealso{
\code{\link{roc.calc}}, \code{\link{combine.roc}}, \code{\link{roc.auc}},\code{\link{roc.plot}}
}
\examples{

true <- c(rep(1, 100), rep(0, 100));
scores <- c( rnorm(100, 1, 1 ), rnorm(100, 0, 1 ) );
roc_mat <- logreg.roc.calc( true, scores );
AUC<- roc.auc(roc_mat);
roc.plot(roc_mat, main=AUC );

}
\keyword{ plot }
