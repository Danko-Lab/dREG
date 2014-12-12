## read_genomic_data -- reads genomic data in the specified zoom parameters...
##

data.one_bigwig<- function(x, chr, bw, depth, window) {
  chr_n <- which(bw$chroms == chr)
#  qbw <- queryByStep.bigWig(bw, bw$chroms[chr_n], 0+x, bw$chromSizes[chr_n], window, do.sum=TRUE)
  qbw <- step.bpQuery.bigWig(bw, bw$chroms[chr_n], 0+x, bw$chromSizes[chr_n], window, with.attributes=FALSE)
  
  if(NROW(qbw) == 0) {
    return(integer(0))
  }
  
  indx <- which(abs(qbw)>depth)
  centers <- (indx-1)*window+x+as.integer(window/2) ## Index is 1-based.  x is the offset.
  return(centers)
}


data.two_bigwig<- function(x, chr, bw1, bw2, depth, window) {
  if(sum(bw1$chroms == chr)==0 | sum(bw2$chroms == chr)==0) {
    ## Add a warning?!
	return(integer(0))
  }

  chrbw1_n <- which(bw1$chroms == chr)
  chrbw2_n <- which(bw2$chroms == chr)

#  indx <- which((queryByStep.bigWig(bw1, bw1$chroms[chrbw1_n], 0+x, bw1$chromSizes[chrbw1_n], window, do.sum=TRUE)>depth) &
#           ((-1)*queryByStep.bigWig(bw2, bw2$chroms[chrbw2_n], 0+x, bw2$chromSizes[chrbw2_n], window, do.sum=TRUE)>depth) ) ## Assume minus values.
  indx <- which((step.bpQuery.bigWig(bw1, bw1$chroms[chrbw1_n], 0+x, bw1$chromSizes[chrbw1_n], window, with.attributes=FALSE)>depth) &
           (step.bpQuery.bigWig(bw2, bw2$chroms[chrbw2_n], 0+x, bw2$chromSizes[chrbw2_n], window, abs.value = TRUE, with.attributes=FALSE)>depth) )


  centers <- (indx-1)*window+x+as.integer(window/2) ## Index is 1-based.  x is the offset.
  return(centers)
}

data.two_bigwig.OR <- function(x, chr, bw1, bw2, depth, window) {
  if(sum(bw1$chroms == chr)==0 & sum(bw2$chroms == chr)==0) { ## This condition is VERY unlikely.
	return(integer(0))
  } else if(sum(bw1$chroms == chr)==0) {
	return(data.one_bigwig(x, chr, bw2, depth, window))
  } else if(sum(bw2$chroms == chr)==0) {
	return(data.one_bigwig(x, chr, bw1, depth, window))
  }

  chrbw1_n <- which(bw1$chroms == chr)
  chrbw2_n <- which(bw2$chroms == chr)
  
#  q1 <- c(queryByStep.bigWig(bw1, bw1$chroms[chrbw1_n], 0+x, bw1$chromSizes[chrbw1_n], window, do.sum=TRUE))
#  q2 <- c(((-1)*queryByStep.bigWig(bw2, bw2$chroms[chrbw2_n], 0+x, bw2$chromSizes[chrbw2_n], window, do.sum=TRUE))) 
  q1 <- c(step.bpQuery.bigWig(bw1, bw1$chroms[chrbw1_n], 0+x, bw1$chromSizes[chrbw1_n], window, with.attributes=FALSE))
  q2 <- c(step.bpQuery.bigWig(bw2, bw2$chroms[chrbw2_n], 0+x, bw2$chromSizes[chrbw2_n], window, abs.value=TRUE, with.attributes=FALSE))

  if(NROW(q1)==0 & NROW(q2)==0) {
    return(integer(0))
  } else if(NROW(q1)==0) {
    indx <- which(q2 > depth)
  } else if(NROW(q2)==0) {
    indx <- which(q1 > depth)
  } else {
    indx <- which(rowSums(data.frame(q1, q2)) >depth) ## Assume minus values.
  }
 
  centers <- (indx-1)*window+x+as.integer(window/2) ## Index is 1-based.  x is the offset.
  return(centers)
}

#' Returns a data.frame with center positions that pass a minimum depth filter ...
#'
#' @param bw Path to bigwig file.
#' @param bw_minus If specified, takes the windows that pass the step in both bigWig files.
#' @param depth Minimum number of reads to return.  
#' @param half_window Distance between to search for #depth reads [bp].
#' @return Returns a data.frame representing a bed file.
get_informative_positions <- function(bw_path, bw_minus_path=NULL, depth= 0, window= 400, step=50, use_OR=TRUE, use_ANDOR=TRUE, debug= TRUE) {
  ## Load bigWigs
  bw  <- load.bigWig(bw_path)
  q_chroms <- bw$chroms[bw$chromSizes > (window+step)]
  
  if(!is.null(bw_minus_path)) {
    bw_minus <- load.bigWig(bw_minus_path)
    q_chroms.minus <- bw_minus$chroms[bw_minus$chromSizes+step > (window+step)]
    q_chroms <- unique(q_chroms, bw_minus$chroms) #chroms[chroms %in% bw_minus$chroms]
  }
  else {
    bw_minus <- NULL
  }

  chroms <- character(0)
  starts <- integer(0)
  for(chr in q_chroms) {
    if(debug) print(chr)
	if(is.null(bw_minus)) {
      vals <- unlist(lapply(seq(0,window,step), data.one_bigwig, chr= chr, bw= bw, depth= depth, window= window))
	}
    else {
	  if(use_ANDOR) {
        windowAND <- 1000; depthAND <- 0; windowOR <- 100; depthOR <- 2
        vals <- c(unlist(lapply(seq(0,window,step), data.two_bigwig.OR, chr= chr, bw1= bw, bw2= bw_minus, depth= depthOR, window= windowOR)), 
                  unlist(lapply(seq(0,window,step), data.two_bigwig,    chr= chr, bw1= bw, bw2= bw_minus, depth= depthAND, window= windowAND)))
	  }
	  else {
        if(use_OR) {
          vals <- unlist(lapply(seq(0,window,step), data.two_bigwig.OR, chr= chr, bw1= bw, bw2= bw_minus, depth= depth, window= window))
        }
        else {
          vals <- unlist(lapply(seq(0,window,step), data.two_bigwig, chr= chr, bw1= bw, bw2= bw_minus, depth= depth, window= window))
        }
      }
    }
	newStarts <- unique(sort(vals))
    starts <- c(starts, newStarts)
	chroms <- c(chroms, rep(chr, NROW(newStarts)))
  }

  ## Unload bigwigs.
  unload.bigWig(bw)
  if(!is.null(bw_minus_path))
    unload.bigWig(bw_minus)
  
  return(data.frame(chrom= chroms, chromStart= starts, chromEnds= starts+1))
}

