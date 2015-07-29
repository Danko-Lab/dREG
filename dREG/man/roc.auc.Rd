\name{roc.auc}
\alias{roc.auc}
\title{ Computes the AUC of a ROC plot.}
\description{
Computes the AUC of a ROC plot.
}
\usage{
roc.auc(ROC)
}
\arguments{
  \item{ROC}{A matrix with 3 columns (FPR, TPR and threshold) calculated by \code{\link{logreg.roc.calc}}. }
}
\details{
The parameter of ROC is a matrix or data frame including 3 columns, FPR(False Positive Rate), TPR(True Positive Rate) and threshold.
}
\value{
AUC value is returned. 
}
\references{
Danko, C. G., Hyland, S. L., Core, L. J., Martins, A. L., Waters, C. T., Lee, H. W., ... & Siepel, A. (2015). 
Identification of active transcriptional regulatory elements from GRO-seq data. Nature methods, 12(5), 433-438.
}
\seealso{
\code{\link{roc.calc}}, \code{\link{logreg.roc.calc}}, \code{\link{combine.roc}}, \code{\link{roc.plot}}
}
\examples{

roc_mat <- data.frame( FPR=c(0, 0.25, 0.5, 0.75, 1), 
	TPR=c(0, 0.5, 0.8, 0.95, 1), 
	threshold=c(1, 1, 1, 1, 1) );
AUC<- roc.auc( roc_mat );
roc.plot( roc_mat, main=AUC );

}
\keyword{ plot }

