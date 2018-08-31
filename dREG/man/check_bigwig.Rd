\name{check_bigwig}
\alias{check_bigwig}
\title{Check bigWig data meet the dREG requirement}
\description{
There are 3 check points for dREG:
1) Each read should be mapped to a locus (3' for PRO-seq and 5' for GRO-seq), not a region;
2) No normalization;
3) Only positive values or negative values in one file;
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
