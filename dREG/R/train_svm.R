## train_svm -- trains an SVM to recognize a certain pattern of regulatory positions.
##

# setClass("regulatory_svm",#"restricted_boltzman_machine",
  # contains="genomic_data_model",
  # representation(
    # asvm= "list"
  # ),
# )


#' Returns a data.frame with center positions that pass a minimum depth filter ...
#'
#' @param gdm Genomic data model.
#' @param bw_plus Path to bigWig file representing the plus strand [char].
#' @param bw_minus Path to bigWig file representing the minus strand [char].
#' @param positions The universe of positions to test and evaluate [data.frame, (chrom,chromCenter)].  Hint: see get_informative_positions().
#' @param positive Bed file containing positive positions [data.frame, (chrom,chromStart,chromEnd)].
#' @param allow Bed file containing positions to avoid in the negative set [data.frame, (chrom,chromStart,chromEnd)].
#' @param n_train Number of training examples.
#' @param n_eval Number of examples on which to test performance.
#' @param pdf_path Specifies the path to a PDF file.  Set to NULL if no PDF should be printed.
#' @param plot_raw_data If TRUE (default), and if a PDF file is specified, plots the raw data used to train the model.
#' @param svm_type "SVR" for support vecctor regression (epsilon-regression).  "P_SVM" for probabilistic SVM (C-classification).
#' @return Returns a trained SVM.
regulatory_svm <- function(gdm, bw_plus_path, bw_minus_path, positions, positive, allow= NULL, n_train=25000, n_eval=1000, pdf_path= "roc_plot.pdf", plot_raw_data=TRUE, extra_enrich_bed= NULL, extra_enrich_frac= 0.1, enrich_negative_near_pos= 0.15, use_rgtsvm=FALSE, svm_type= "SVR", ..., debug= TRUE) {
  ########################################
  ## Divide into positives and negatives.

  if (use_rgtsvm)
  {
    if(!requireNamespace("Rgtsvm"))
      stop("Rgtsvm has not been installed fotr GPU computing.");

    predict = Rgtsvm::predict.gtsvm;
    svm = Rgtsvm::svm;
  }

  if( class(asvm)=="svm" && use_rgtsvm) class(asvm)<-"gtsvm";
  if( class(asvm)=="gtsvm" && !use_rgtsvm) class(asvm)<-"svm";

  inter_indx <- (n_train+n_eval)
  indx_train <- c(1:n_train, (inter_indx+1):(inter_indx+n_train))
  indx_eval  <- c((n_train+1):(inter_indx), (inter_indx+n_train+1):(2*inter_indx))

  ## Read genomic data.
  if(debug) print("Collecting training data.")
  if(length(bw_plus_path) == 1) {
    tset <- get_test_set(positions= positions, positive= positive, allow= allow, n_samp= (n_train+n_eval), extra_enrich_bed= extra_enrich_bed, extra_enrich_frac= extra_enrich_frac, enrich_negative_near_pos= enrich_negative_near_pos)

    ## Get training indices.
    x_train_bed <- tset[indx_train,c(1:3)]
    y_train <- tset[indx_train,4]
    x_predict_bed <- tset[indx_eval,c(1:3)]
    y_predict <- tset[indx_eval,4]

	## Write out a bed of training positions to avoid during test ...
    if(debug) {
      write.table(x_train_bed, "TrainingSet.bed", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
	  write.table(indx_train, "TrainIndx.Rflat")
    }

    x_train <- read_genomic_data(gdm, x_train_bed, bw_plus_path, bw_minus_path)
  } else {
    x_train <- NULL
    y_train <- NULL
    stopifnot(NROW(bw_plus_path) == NROW(bw_minus_path) & NROW(bw_plus_path) == NROW(positive))
    for(x in 1:length(bw_plus_path)){
      tset_x <- get_test_set(positions= positions[[x]], positive= positive[[x]], allow= allow[[x]], n_samp= (n_train+n_eval), extra_enrich_bed= extra_enrich_bed[[x]], extra_enrich_frac= extra_enrich_frac, enrich_negative_near_pos= enrich_negative_near_pos)

      x_train_bed <- tset_x[indx_train,c(1:3)]
      y_train <- c(y_train, tset_x[indx_train,4])

      x_train <- rbind(x_train, read_genomic_data(gdm, x_train_bed, bw_plus_path[[x]], bw_minus_path[[x]]))
    }
  }

  ########################################
  ## Train the model.
  if(debug) print("Fitting SVM.")
  if (svm_type == "SVR") {
    if(debug) print("Training a epsilon-regression SVR.")
    asvm <- svm( x_train, y_train )
  }
  if (svm_type == "P_SVM") {
    if(debug) print("Training a probabilistic SVM.")
    asvm <- svm( x_train, as.factor(y_train), probability=TRUE)
  }

  ########################################
  ## If a PDF file is specified, test performance, and write ROC plots to a PDF file.
  ## Currently *NOT* supported when training with >1 dataset.
 if(!is.null(pdf_path) && !is.na(pdf_path) && n_eval>0 && length(bw_plus_path) == 1) {
  pdf(pdf_path)
    # Plot raw data, if desired.
    if(plot_raw_data) {
      plot(colSums(x_train[y_train == 1,]), ylab="Training data", type="l", ...)
      points(colSums(x_train[y_train == 0,]), col="gray", type="l", ...)
    }
    remove(x_train)

    ## Predict on a randomly chosen set of sequences.
    if(debug) print("Collecting predicted data.")

    x_predict <- read_genomic_data(gdm, x_predict_bed, bw_plus_path, bw_minus_path)

    pred <- predict( asvm, x_predict )

    ## Plot raw prediction data, if desired.
    if(plot_raw_data) {
      plot(colSums(x_predict[y_predict == 1,]), ylab="Prediction data", type="l", ...)
      points(colSums(x_predict[y_predict == 0,]), col="gray", type="l", ...)
    }
    remove(x_predict)

    ## Write ROC plots.
    roc_values <- logreg.roc.calc(y_predict, pred)
    AUC<- roc.auc(roc_values)
    roc.plot(roc_values, main=AUC, ...)
    print(paste("Model AUC: ",AUC))
    remove(roc_values)
  dev.off()
 }
 else {
   remove(x_train)
 }

  return(asvm)
}


