#' Returns a list of two data.frames first representing the _train set, second representing the _eval set ...
#'
#' @param positions The universe of positions to test and evaluate [data.frame, (chrom,chromCenter)].  Hint: see get_informative_positions().
#' @param positive Bed file containing positive positions [data.frame, (chrom,chromStart,chromEnd)].
#' @param n_train Number of training examples.
#' @param allow    Bed file containing inverse negative set of positions [data.frame, (chrom,chromStart,chromEnd)].
#' @param enrich_negative_near_pos Fraction of training examples chosen to be nearby (<=5kb) a positive example [0,1].  Default= 0.25.
#' @param extra_enrich_bed Extra bed file to enrich near (default=NULL).
#' @param extra_enrich_frac Fraction of final positions sampled in the negative set which are in the bed file.  Unused if extra_enrich_bed is NULL.
#' @param n_eval Number of examples on which to test performance.
#' @return Returns a list of two data.frames first representing the _train set, second representing the _eval set.

get_test_set <- function(positions, positive, n_samp, allow=NULL, enrich_negative_near_pos= 0.15, extra_enrich_bed= NULL, extra_enrich_frac= 0.1, avoid_dist= 100) {
  if(enrich_negative_near_pos < 0 | enrich_negative_near_pos > 1) stop('ERROR: enrich_negative_near_pos must be in the rage [0,1]!')
  if(is.null(extra_enrich_bed)) {
    extra_enrich_frac=0
  }
  else {
    if(extra_enrich_frac < 0 | extra_enrich_frac > 1) stop('ERROR: extra_enrich_frac must be in the rage [0,1]!')
  }
  
  ########################################
  ## Divide into positives and negatives.
  all_feat  <- feat(seqname= positions[,1], start= positions[,2], end= (positions[,3]))
  positive_feat <- feat(seqname= positive[,1], start= positive[,2], end= positive[,3])
  ol <- overlap.feat(x= all_feat, filter= positive_feat)
  pos_indx <- match(paste(ol$seqname, ol$start, ol$end), paste(all_feat$seqname, all_feat$start, all_feat$end))
  
  ## Get negatives.  Make sure negatives are not right beside positives.  Choice of 100 bp is arbitrary ... longer windows might widen regions; shorter might make training more challenging, and increase the apparent FPR.
  avoid <- positive
  
## NOTE: This line will crap out of the column names are different.  
  if(!is.null(allow)) avoid <- rbind(avoid[,c(1:3)], allow[,c(1:3)]) 
  
  avoid[,2] <- avoid[,2]-avoid_dist
  avoid[,3] <- avoid[,3]+avoid_dist
  avoid_feat <- feat(seqname= avoid[,1], start= avoid[,2], end= avoid[,3])
  ol <- overlap.feat(x= all_feat, filter= avoid_feat)
  avoid_indx <- match(paste(ol$seqname, ol$start, ol$end), paste(all_feat$seqname, all_feat$start, all_feat$end))

  neg_indx <- rep(TRUE, NROW(all_feat))
  neg_indx[avoid_indx] <- FALSE
  neg_indx <- which(neg_indx)
  n_samp_neg <- (1-enrich_negative_near_pos-extra_enrich_frac)*n_samp
  n_samp_extra_enrich <- extra_enrich_frac*n_samp
  
  ## Subsample.
  pos_indx_train_test <- sample(pos_indx, n_samp, replace=FALSE)
  neg_indx_train_test <- sample(neg_indx, n_samp_neg, replace=FALSE)
    
  ## ID regions near positives.
  if(enrich_negative_near_pos > 0) {
   near_size <- 5000 #[bp]
   nearby_feat <- feat(seqname=positive[,1], start= (positive[,2]-near_size), end= (positive[,3]+near_size))
   ol <- overlap.feat(x= all_feat, filter= nearby_feat)
   nearby_indx <- match(paste(ol$seqname, ol$start, ol$end), paste(all_feat$seqname, all_feat$start, all_feat$end))
   nearby_bool <- rep(FALSE, NROW(all_feat)) ## Start FALSE
   nearby_bool[nearby_indx] <- TRUE  ## Turn just those nearby to TRUE
   nearby_bool[avoid_indx] <- FALSE ## Turn those inside postive regions to FALSE
   nearby_bool[neg_indx_train_test] <- FALSE ## Turn those which have already been selected to FALSE
   nearby_bool <- which(nearby_bool)
   n_samp_nearby <- n_samp-n_samp_neg-n_samp_extra_enrich
   neg_indx_train_test <- sort(c(neg_indx_train_test, sample(nearby_bool, n_samp_nearby, replace=FALSE)))
  }
  
  ## ID regions inside of an extra 'enrich' bed file set.
  if(!is.null(extra_enrich_bed)) {
   extra_feat <- feat(seqname=extra_enrich_bed[,1], start= extra_enrich_bed[,2], end= extra_enrich_bed[,3])
   ol <- overlap.feat(x= all_feat, filter= extra_feat)
   extra_indx <- match(paste(ol$seqname, ol$start, ol$end), paste(all_feat$seqname, all_feat$start, all_feat$end))
   extra_bool <- rep(FALSE, NROW(all_feat)) ## Start FALSE
   extra_bool[extra_indx] <- TRUE  ## Turn just those in extra regions to TRUE
   extra_bool[avoid_indx] <- FALSE ## Turn those inside postive regions to FALSE
   extra_bool[neg_indx_train_test] <- FALSE ## Turn those which have already been selected to FALSE   
   extra_bool <- which(extra_bool)
   neg_indx_train_test <- sort(c(neg_indx_train_test, sample(extra_bool, n_samp_extra_enrich, replace=FALSE)))
  }

  return(cbind(
	positions[c(pos_indx_train_test[c(1:n_samp)], neg_indx_train_test[c(1:n_samp)]),], 
	c(rep(1,n_samp), rep(0,n_samp)))) ## Train on the first n_train examples.
}