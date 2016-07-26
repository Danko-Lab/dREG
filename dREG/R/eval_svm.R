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
eval_reg_svm <- function(gdm, asvm, positions, bw_plus_path, bw_minus_path, batch_size=50000, ncores=3, use_rgtsvm=FALSE, debug= TRUE) {

  if(batch_size>NROW(positions)) batch_size= NROW(positions)
  n_elem <- NROW(positions)
  n_batches <- floor(n_elem/batch_size)
  interval <- unique(c( seq( 1, n_elem+1, by = batch_size ), n_elem+1))

  pos.ord <- order(positions[,1], positions[,2], positions[,3]);
  pos.sorted <- positions[pos.ord,];

  do.extract<-function(x)
  {
     batch_indx<- c( interval[x]:(interval[x+1]-1) );
     return( read_genomic_data(gdm, pos.sorted[batch_indx,,drop=F], bw_plus_path, bw_minus_path) );
  }

  do.predict <- function( mat_features ){
     if(asvm$type == 0) { ## Probabilistic SVM
       batch_pred <- svm_predict(asvm, mat_features, use_rgtsvm=use_rgtsvm, probability=TRUE)
     }
     else { ## epsilon-regression (SVR)
       batch_pred <- svm_predict(asvm,mat_features, use_rgtsvm=use_rgtsvm, probability=FALSE)
     }
     return(batch_pred);
  }
 
  ## Do elements of each intervals
  if(!use_rgtsvm)
  {
    scores<- unlist(mclapply(c(1:(length(interval)-1)), function(x) {
      print(paste(x, "of", length(interval)-1) );
      pred <- do.predict( do.extract(x) );
      gc(); 
      return( pred );
    }, mc.cores= ncores))
  }else
  {
    n.loop <- ceiling((length(interval)-1)/ncores);
     scores <- unlist( lapply(1:n.loop, function(i) {
       n.start = (i-1)*ncores+1;
       n.stop = ifelse( length(interval)-1 <= i*ncores, length(interval)-1, i*ncores );
       feature_list<- mclapply(n.start:n.stop, function(x) {
          print(paste(x, "of", length(interval)-1) );
          return( do.extract(x) );
       }, mc.cores= ncores);
      
      pred <- do.predict( do.call("rbind", feature_list) );
      feature_list <- NULL;
      gc();
      return( pred );
     } ));
  }

  ## Test code
  ## all( pos.sorted[ order(pos.ord),  ] == positions );

  ## sort back the genome loci.
  scores <- scores[ order(pos.ord) ];

  return(as.double(scores))
}

svm_predict<-function( asvm, mat_features, use_rgtsvm=FALSE, probability=FALSE )
{
	if (use_rgtsvm)
		require(Rgtsvm);

	if( class(asvm)=="svm" && use_rgtsvm) class(asvm)<-"gtsvm";
	if( class(asvm)=="gtsvm" && !use_rgtsvm) class(asvm)<-"svm";

	ret <- c();
	if(use_rgtsvm)
		ret <- Rgtsvm::predict.gtsvm(asvm, mat_features, probability=probability)
	else
		ret <- predict( asvm, mat_features, probability=probability );

	return(ret);
}

