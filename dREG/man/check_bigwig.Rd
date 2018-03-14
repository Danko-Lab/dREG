\name{check_bigwig}
\alias{check_bigwig}
\title{Check bigWig data meet the dREG requirement}
\description{
There are 3 check points for dREG:
1) No normalization;
2) Positive values in plus strand and negative values in minus strand ;
3) Each read should be mapped to a locus, not a region;
}
\usage{
check_bigwig(bw_path, strand = "+", out.file = "")
}
\arguments{
  \item{bw_path}{String value, bigWig file}
  \item{strand}{"+" or "-", strand}
  \item{out.file}{file name, indicating detailed information will be outputted.}
}
\value{
Boolean value indicates whether it is suitable to do peak calling. If the bigWig doesn't meet the requirements of dREG, the function will return FALSE with the details outputted into console or file. 
}
