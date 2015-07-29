\name{combine.roc}
\alias{combine.roc}
\title{ Combines ROC plots}
\description{Combines ROC plots, interpolating and weighting by nTP.}
\usage{
combine.roc(list.roc, 
      weight = rep(1, NROW(list.roc)), 
      interp.corners = FALSE, 
      use.max = FALSE, 
      nvals = 100)
}

\arguments{
  \item{list.roc}{List including multiple ROC data frame}
  \item{weight}{Weight vector for each ROC dataframe}
  \item{interp.corners}{Logical value indicating if the header(1,1) and tail values(0,0) are interpolated to each ROC data frame.}
  \item{use.max}{Logical value indicating if maximum value of muliple ROCs at same point are used as TPF values.}
  \item{nvals}{Integer value indicating interval number for ROC plot.}
}
\details{
}
\value{
A data frame with 2 columns is returned
\item{FPR}{False Positive Rate}
\item{TPR}{True Positive Rate}
}
\references{
Danko, C. G., Hyland, S. L., Core, L. J., Martins, A. L., Waters, C. T., Lee, H. W., ... & Siepel, A. (2015). 
Identification of active transcriptional regulatory elements from GRO-seq data. Nature methods, 12(5), 433-438.
}

\seealso{
\code{\link{roc.calc}}, \code{\link{logreg.roc.calc}}, \code{\link{roc.auc}}, \code{\link{roc.plot}}
}

\examples{

list.roc<-list();

true <- c(rep(1, 100), rep(0, 100));
scores <- c( rnorm(100, 1, 1 ), rnorm(100, 0, 1 ) );
list.roc[[1]] <- logreg.roc.calc( true, scores );

true <- c(rep(1, 120), rep(0, 110));
scores <- c( rnorm(120, 1, 0.8 ), rnorm(110, 0, 1.2 ) );
list.roc[[2]] <- logreg.roc.calc( true, scores );

r <- combine.roc(list.roc);
roc.plot(r)

}
\keyword{ plot }
