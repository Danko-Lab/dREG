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
      use_rgtsvm = TRUE, 
      min_score = NULL, 
      pv_adjust="fdr",
      pv_threshold=0.05,
      smoothwidth = 4,
      cpu_cores=1, 
      gpu_cores=1)
}
\arguments{
  \item{asvm}{SVR model pre-trained for dREG package, which can be downloaded from the dREG page in Github.}
  \item{gdm}{Genomic data model return by \code{\link{genomic_data_model}}. This data is binding with SVR model.}
  \item{bw_plus_path}{File name indicating file path to bigWig file representing the plus strand.}
  \item{bw_minus_path}{File name  indicating file path to bigWig file representing the minus strand.}
  \item{infp_bed}{A BED data frame indicating informative sites and scores returned by \code{\link{eval_reg_svm}}. If NULL is specified, the peak calling starts from the informative sites finding and predicting. }
  \item{use_rgtsvm}{Logical value indictating whether the predict will be performed on GPU through the Rgtsvm package.}
  \item{min_score}{Numerical value indicating the minimum dREG score applied to the peak calling procedure. If NULL is specified, this value is calculated based on the predicted scores. }
  \item{pv_adjust}{String value indictating which correction method is used to do multiple comparison, see details in \code{\link{p.adjust}}, default is 'fdr'. }
  \item{pv_threshold}{Numerical value indicating the threshold is used to report the dREG peaks.}
  \item{smoothwidth}{Numerical value indicating the parameter of curve smooth in the moving average.}
  \item{cpu_cores}{Number of CPU cores in parallel computing.}
  \item{gpu_cores}{Number of GPU cores in parallel computing if \pkg{Rgtsvm} is used.}
}
\details{

}
\value{

This function returns a list containing 6 items, including: \cr

1) dREG peaks: \code{peak_bed} \cr

\item{chr}{ Chromosome}
\item{start}{ Start position}
\item{end}{End position}
\item{score}{Maxmimum score in the peak region}
\item{prob}{Probability of multivariate Laplace distribution indicating the probability of the peak points belonging to negative set (No divergent peak). }
\item{center}{the center position in original peak}

2) Informative sites with score \code{infp_bed} \cr
\item{chr}{ Chromosome}
\item{start}{ Start position}
\item{end}{End position}
\item{score}{ predicted score}
\item{infp}{ indicating the informative site or dense site}

3) Broad peak regions \code{peak_broad}
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

5) Raw results of peak calling \code{raw_peak}


6) Summary of peak calling \code{peak_sum}
\item{adjust.none.0.05}{The count of dREG peak withou p-value correction.}
\item{adjust.fdr.0.05}{The count of dREG peak adjusted by the 'fdr' method.}
\item{adjust.BH.0.05}{The count of dREG peak adjusted by the 'BH' method.}
\item{adjust.bonferroni.0.05}{The count of dREG peak adjusted by the 'bonferroni' method.}
\item{adjust.holm.0.05}{The count of dREG peak adjusted by the 'holm' method.}
\item{adjust.hochberg.0.05}{The count of dREG peak adjusted by the 'hochberg' method.}
\item{adjust.BY.0.05}{The count of dREG peak adjusted by the 'BY' method.}
\item{peak.sig.score}{The score range of significant dREG peaks.}
\item{peak.narrow100}{The ignored narrow peaks which length are less than 100.}
\item{peak.narrow100.sig}{The ignored narrow peaks which may be significant based on 'peak.sig.score'.}
\item{peak.narrow100.score}{The score range of the ignored narrow peaks.}

}
\examples{

# load("../asvm.6.6M.20170828.rdata");
# gdm <- genomic_data_model(window_sizes= c(10, 25, 50, 500, 5000), half_nWindows= c(10, 10, 30, 20, 20) )
# bw_plus_path <- "K562.chr21.plus.bw"
# bw_minus_path <- "K562.chr21.minus.bw"
# x <- peak_calling( svm, gdm, bw_plus_path, bw_minus_path, cpu_cores=12, use_rgtsvm=T)
# show(x$peak_bed);


}
\keyword{ function }

