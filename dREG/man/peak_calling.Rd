\name{peak_calling}
\alias{peak_calling}
\title{
Peak calling based on dREG prediction
}
\description{
This procedure calls SVR prediction for paired bigWig files using pre-trained SVM model and detects divergent peaks based on the predicted score.
}
\usage{
peak_calling(asvm, gdm, 
      bw_plus_path, bw_minus_path, 
      infp_bed = NULL, 
      ncores = 1, 
      use_rgtsvm = TRUE, 
      min_score = NULL, 
      smoothwidth = 4)
}
\arguments{
  \item{asvm}{SVR model pre-trained for dREG package, which can be downloaded from the dREG page in Github.}
  \item{gdm}{Genomic data model return by \code{\link{genomic_data_model}}. This data is binding with SVR model.}
  \item{bw_plus_path}{File name indicating file path to bigWig file representing the plus strand.}
  \item{bw_minus_path}{File name  indicating file path to bigWig file representing the minus strand.}
  \item{infp_bed}{A BED data frame indicating informative sites and scores returned by \code{\link{eval_reg_svm}}. If NULL is specified, the peak calling starts from the informative sites finding and predicting. }
  \item{ncores}{Number of CPU cores in parallel computing.}
  \item{use_rgtsvm}{Logical value indictating whether the predict will be performed on GPU through the Rgtsvm package.}
  \item{min_score}{Numerical value indicating the minimum dREG score applied to the peak calling procedure. If NULL is specified, this value is calculated based on the predicted scores. }
  \item{smoothwidth}{Numerical value indicating the parameter of curve smooth in the moving average.}
}
\details{

}
\value{

This function returns a list containing 3 items, which include: \cr

1) dREG Peaks: \code{peak_bed} \cr

\item{chr}{ Chromosome}
\item{start}{ Start position}
\item{end}{End position}
\item{score}{Maxmimum score in the peak region}
\item{prob}{Probability of multivariate Laplace distribution indicating the probability of the peak points belonging to negative set (No divergent peak). }
\item{center}{the center position in original peak}

2) Informative Sites with score \code{infp_bed} \cr
\item{chr}{ Chromosome}
\item{start}{ Start position}
\item{end}{End position}
\item{score}{ predicted score}
\item{infp}{ indicating the informative site or dense site}

3) Summary of peak region \code{peak_sum}
\item{chr}{ Chromosome}
\item{start}{ Start position}
\item{end}{End position}
\item{no}{ index }
\item{min}{ minimum score in this region}
\item{max}{ maximum score in this region}
\item{mean}{ score mean in this region}
\item{sum}{ score um in this region }
\item{stdev}{ standard deviation of scores in this region}
\item{count}{ informative site in this region}

4) Threshold of dREG score \code{min_score}



}
\examples{

# load("../asvm.6.6M.20170828.rdata");
# gdm <- genomic_data_model(window_sizes= c(10, 25, 50, 500, 5000), half_nWindows= c(10, 10, 30, 20, 20) )
# bw_plus_path <- "K562.chr21.plus.bw"
# bw_minus_path <- "K562.chr21.minus.bw"
# x <- peak_calling( svm, gdm, bw_plus_path, bw_minus_path, ncores=12, use_rgtsvm=T)
# show(x$peak_bed);


}
\keyword{ function }

