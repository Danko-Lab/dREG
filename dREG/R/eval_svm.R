## train_svm -- trains an SVM to recognize a certain pattern of regulatory positions.
##

#' Evaluates a set of genomic coordinates for regulatory potential using P/GRO-seq data ...
#'
#' @param asvm A pre-trained SVM model from the e1071 package.
#' @param positions The universe of positions to test and evaluate [data.frame, (chrom,chromCenter)].  Hint: see get_informative_positions().
#' @param bw_plus_path Path to bigWig file representing the plus strand [char].
#' @param bw_minus_path Path to bigWig file representing the minus strand [char].
#' @param batch_size Number of positions to evaluate at once (more might be faster, but takes more memory).
#' @return Returns the value of the SVM for each genomic coordinate specified.
eval_reg_svm <- function(gdm, asvm, positions, bw_plus_path, bw_minus_path, batch_size=50000, ncores=3, debug= TRUE) {

  if(batch_size>NROW(positions)) batch_size= NROW(positions)
  n_elem <- NROW(positions)
  n_batches <- floor(n_elem/batch_size)
  interval <- unique(c( seq( 1, n_elem+1, by = batch_size ), n_elem+1))

  ## Do elements of each intervals
  scores<- unlist(mclapply(c(1:(length(interval)-1)), function(x) {

    print(paste(x, "of", n_batches))
    batch_indx<- c( interval[x]:(interval[i+1]-1) )
    x_predict <- read_genomic_data(gdm, positions[batch_indx,,drop=F], bw_plus_path, bw_minus_path)
    if(asvm$type == 0) { ## Probabilistic SVM
      batch_pred <- predict(asvm, x_predict, probability=TRUE)
    }
    else { ## epsilon-regression (SVR)
      batch_pred <- predict(asvm, x_predict)
    }
    return(batch_pred)
  }, mc.cores= ncores))

  return(as.double(scores))
}
