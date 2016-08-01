## read_genomic_data -- reads genomic data in the specified zoom parameters...
##

setClass("genomic_data_model",#"restricted_boltzman_machine",
  representation(
    n_zooms="integer",
    window_sizes="integer",
    half_nWindows="integer"
    ),
)

genomic_data_model <- function(window_sizes, half_nWindows) {
  stopifnot(NROW(window_sizes) == NROW(half_nWindows))
  new("genomic_data_model", n_zooms= as.integer(NROW(window_sizes)), window_sizes= as.integer(window_sizes), half_nWindows= as.integer(half_nWindows))
}

#' Reads genomic data from the specified position...
#'
#' @param bed A data.frame of genomic regions in bed format (chrom, start, end).
#' @param bigwig_plus Path to bigwig file representing GRO-seq/ PRO-seq reads on the plus strand.
#' @param bigwig_minus Path to bigwig file representing GRO-seq/ PRO-seq reads on the minus strand.
#' @param as_matrix If true, returns a matrix object.
#' @param ncores The number of cores.
#' @param scale.method Default is logistic, but if set to linear it will return read counts normalized by total read count
#' @return Returns a list() object, where each element in the list is the zoom data
#' centered on a
read_genomic_data <- function(gdm, bed, file_bigwig_plus, file_bigwig_minus, as_matrix= TRUE, scale.method=c("logistic", "linear"), batch_size=50000, ncores=1) {

  stopifnot(NROW(gdm@window_sizes) == NROW(gdm@half_nWindows))

  zoom<- list(as.integer(gdm@window_sizes), as.integer(gdm@half_nWindows))
  #batch_size = 50000;
  n_elem = NROW(bed)
  n_batches = floor( n_elem/batch_size )
  if(n_batches < ncores)
  {
	  batch_size = ceiling( n_elem/ncores );
	  n_batches = floor( n_elem/batch_size);
  }

  interval <- unique(c( seq( 1, n_elem+1, by = batch_size ), n_elem+1))

  if(missing(scale.method)){ scale.method <- "logistic" };

  total.read.count<- sum(abs(get_reads_from_bigwig(file_bigwig_plus, file_bigwig_minus)));

  bed.ord <- order(bed[,1], bed[,2], bed[,3]);
  if ( all( bed.ord == c(1:NROW(bed)) ) )
     bed.sorted <- bed
  else
     bed.sorted <- bed[ bed.ord, ];

  datList <- list();
  for(i in 1:ceiling( (length(interval)-1)/ncores ))
  {
	  start_batch <- (i-1)*ncores+1;
	  stop_batch <- i*ncores;
	  if(stop_batch>(length(interval)-1)) stop_batch <- length(interval)-1;

	  datList[start_batch:stop_batch] <- mclapply(start_batch:stop_batch, function(x) {
		  batch_indx<- c( interval[x]:(interval[x+1]-1) )

		  # The output from C/C++ is changed to matrix(n_windows, n_sample)) since 6/20/2016
		  # The original result was a list, which needs to rbind() to matrix.

		  if(scale.method=="logistic"){
			  dat <- .Call("get_genomic_data_R",
							as.character( bed.sorted[ batch_indx,1 ] ),
							as.integer( floor((bed.sorted[ batch_indx,3 ] + bed.sorted[ batch_indx,2 ])/2) ),
							as.character( file_bigwig_plus ),
							as.character( file_bigwig_minus ),
							zoom,
							as.logical(TRUE),
							PACKAGE= "dREG")
		  }
		  else{
			  dat <- .Call("get_genomic_data_R",
						   as.character( bed.sorted[ batch_indx,1 ] ),
						   as.integer( floor((bed.sorted[ batch_indx,3 ] + bed.sorted[ batch_indx,2 ])/2) ),
						   as.character( file_bigwig_plus ),
						   as.character( file_bigwig_minus ),
						   zoom,
						   as.logical(FALSE),
						   PACKAGE= "dREG")

			  if( !is.null(dat) )
				  ##dat<-lapply(dat, "/", total.read.count);
				  dat <- dat/total.read.count;
		  }

		  if( is.null(dat))
		  	stop("Failed to Call C/C++ functions.\n");

#	cat(x, NROW(dat), NCOL(dat), "\n");
		  return( as.data.frame(t(dat) ) );
	  }, mc.cores=ncores);
  }

  if( length(datList)==1)
    dat <- datList[[1]]
  else
  {
    if(requireNamespace("data.table"))
	  dat <- as.matrix( data.table::rbindlist(datList) )
    else
  	  dat <- as.matrix( do.call(rbind, datList) );
  }

  rm(datList);

  if ( !all( bed.ord == c(1:NROW(bed)) ) )
     ## all(bed.sorted[order(bed.ord),] == bed)
     dat <- dat [ order(bed.ord), ];

  if( !as_matrix )
	 dat <- c(t(dat));

  return(dat);
}

# query read counts of all chromosomes from bigWig file.
#
# @bw.plus,     bigWig filename
# @bw.minus,    bigWig filename
# @chromInfo    data.frame with 2 columns(chr, size);
#
# @return       vector of reads in plus and minus file.
get_reads_from_bigwig <- function(file_bigwig_plus, file_bigwig_minus)
{
    bw.plus  <- load.bigWig(file_bigwig_plus)
    bw.minus <- load.bigWig(file_bigwig_minus)

	## 1) It takes long time
	## 2) The offset is too big for some unmapped section, it will cause errors in library bigWig

    #offset_dist <- 250;
    #df.bed.plus<-data.frame(bw.plus$chroms, offset_dist, bw.plus$chromSizes, names=".", scores=".",strands="+")
    #df.bed.minus<-data.frame(bw.minus$chroms, offset_dist, bw.minus$chromSizes, names=".", scores=".", strands="-")
    #r.plus <- sum(abs(bed6.region.bpQuery.bigWig( bw.plus, bw.minus, df.bed.plus)));
    #r.minus <- sum(abs(bed6.region.bpQuery.bigWig( bw.plus, bw.minus, df.bed.minus)));

 	r.plus  <- round(bw.plus$mean * bw.plus$basesCovered );
 	r.minus <- round(bw.minus$mean * bw.minus$basesCovered );

    try( unload.bigWig( bw.plus ) );
    try( unload.bigWig( bw.minus ) );

    return(c(r.plus,r.minus));
}