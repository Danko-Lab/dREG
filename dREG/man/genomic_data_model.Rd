\name{genomic_data_model}
\alias{genomic_data_model}
\title{
Builds a multiscale window object for feature extracting.
}
\description{
}
\usage{
genomic_data_model(window_sizes, half_nWindows)
}
\arguments{
  \item{window_sizes}{vector of integer, indicating the genomic size (bp) for each window.}
  \item{half_nWindows}{vector of integer, specifying the window count for each above window. Because the windows are extended at the both sides of an observed position, here this number is considered as half number(left or right side). }
}
\details{
The total number of features including plus and strand is sum(half_nWindows)*2 sides * 2 strands.\cr
The covered region are max(window_sizes)*half_nWindows[which.max(window_sizes)]*2 bps.
}
\value{
A S4 object including two attributes.
}
\examples{
  gdm <- genomic_data_model(window_sizes = c(10, 25, 50, 500, 5000), 
                            half_nWindows= c(10, 10, 30, 20, 20) );
}
