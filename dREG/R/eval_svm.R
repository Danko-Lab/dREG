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

  if(!file.exists(bw_plus_path))
    stop( paste("Can't find the bigwig of plus strand(", bw_plus_path, ")"));

  if(!file.exists(bw_minus_path))
    stop( paste("Can't find the bigwig of minus strand(", bw_minus_path, ")"));

  if( batch_size>NROW(positions))
      batch_size= NROW(positions)

  if(NROW(positions)/ncores < batch_size)
      batch_size <- ceiling(NROW(positions)/ncores);

  n_elem <- NROW(positions)
  interval <- unique(c( seq( 1, n_elem+1, by = batch_size ), n_elem+1))

  pos.ord <- order(positions[,1], positions[,2], positions[,3]);
  pos.sorted <- positions[pos.ord,];

  do.predict <- function( mat_features){

    if (use_rgtsvm)
    {
      if(!requireNamespace("Rgtsvm"))
        stop("Rgtsvm has not been installed fotr GPU computing.");

      predict = Rgtsvm::predict.gtsvm;
    }

    if( class(asvm)=="svm" && use_rgtsvm) class(asvm)<-"gtsvm";
    if( class(asvm)=="gtsvm" && !use_rgtsvm) class(asvm)<-"svm";

    ## if using the preload way to accerate the predict speed
    if( !is.null(asvm$cluster) || !is.null(asvm$pointer) )
    {
        if(asvm$type == 0) { ## Probabilistic SVM
          batch_pred <- Rgtsvm::predict.run( asvm, mat_features, probability=TRUE );
        }
        else { ## epsilon-regression (SVR)
          batch_pred <- Rgtsvm::predict.run( asvm, mat_features, probability=FALSE)
        }
    }
    else
    {
        if(asvm$type == 0) { ## Probabilistic SVM
          batch_pred <- predict( asvm, mat_features, probability=TRUE );
        }
        else { ## epsilon-regression (SVR)
          batch_pred <- predict( asvm, mat_features, probability=FALSE)
        }
    }

    return(batch_pred);
  }

  ## Do elements of each intervals
  if(!use_rgtsvm)
  {
     cpu.fun <- function(x)
     {
        require("dREG");
          batch_idx <- c( interval[x]:(interval[x+1]-1) );
          feature <- read_genomic_data(gdm, pos.sorted[batch_idx,,drop=F], bw_plus_path, bw_minus_path)
          pred <- do.predict( feature );
          gc();
          return( pred );
     }

    sfInit(parallel = TRUE, cpus = ncores, type = "SOCK" )
    sfExport("gdm", "pos.sorted", "bw_plus_path", "bw_minus_path", "interval", "asvm", "do.predict", "use_rgtsvm");

    fun <- as.function(cpu.fun);
    environment(fun)<-globalenv();

    scores <- unlist( sfLapply( 1:(length(interval)-1), fun) );
    sfStop();

  }
  else
  {
    n.loop <- ceiling((length(interval)-1)/ncores);

    scores <- unlist( lapply(1:n.loop, function(i) {
      n.start = (i-1)*ncores+1;
      n.stop = ifelse( length(interval)-1 <= i*ncores, length(interval)-1, i*ncores );

      #if(!use_snowfall)
      #{
      #   feature_list<- mclapply(n.start:n.stop, function(x) {
      #      print(paste(x, "of", length(interval)-1) );
      #      batch_indx<- c( interval[x]:(interval[x+1]-1) );
      #      return(read_genomic_data(gdm, pos.sorted[batch_indx,,drop=F], bw_plus_path, bw_minus_path));
      #   }, mc.cores= ncores);
      #}
      #else
      #{
      if(!requireNamespace("snowfall"))
         stop("Snowfall has not been installed fotr big data.");

      pos.list = list();
      for(x in n.start:n.stop)
      {
        print(paste(x, "of", length(interval)-1) );
        batch_indx<- c( interval[x]:(interval[x+1]-1) );
        pos.list[[x-n.start+1]] <- pos.sorted[batch_indx,,drop=F];
      }

      cpu.fun <- function(pos.bed)
      {
        require("dREG");

        cat("PID=", Sys.getpid(), "\n");
        return(read_genomic_data(gdm, pos.bed, bw_plus_path, bw_minus_path));
      }

      sfInit(parallel = TRUE, cpus = ncores, type = "SOCK" )
      sfExport("gdm", "bw_plus_path", "bw_minus_path");

      fun <- as.function(cpu.fun);
      environment(fun)<-globalenv();

      feature_list <- list();
      if(length(pos.list)>1)
         feature_list <- sfLapply( pos.list, fun)
      else
         feature_list[[1]] <- cpu.fun(pos.list[[1]]);
      sfStop();
     #}

      pred <- do.predict( do.call("rbind", feature_list) );

      rm( feature_list );
      gc(verbose=T, reset=T);

      return( pred );
     } ));
  }

  ## Test code
  ## all( pos.sorted[ order(pos.ord),  ] == positions );

  ## sort back the genome loci.
  scores <- scores[ order(pos.ord) ];

  return(as.double(scores))
}

