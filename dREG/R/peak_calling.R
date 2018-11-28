peak_calling_nopred<-function( infp_bed, min_score, pv_adjust="fdr", pv_threshold=0.05, smoothwidth=4, ncores=1 )
{
  colnames(infp_bed) <- c("chr", "start", "end", "score", "infp");

  broadpeak_sum <- get_broadpeak_summary( infp_bed[,-5], threshold=0.05 );
  if (NROW(broadpeak_sum) != sum( !is.na(broadpeak_sum$max)) )
  {
    ## remove unknown contig, e.g.g chr1_gl000192_random chr10, chr17_ctg5_hap1, chr6_mann_hap4,chr6_qbl_hap6
    ## Notice: chr1_gl000192_random causes Bedmap is failed to get summary information for chr10-chr19
    ##         remove these contigs temporally
    infp_bed <- infp_bed[ grep("_", infp_bed$chr,invert=TRUE), ]
    broadpeak_sum <- get_broadpeak_summary( infp_bed[,-5], threshold=0.05 );
  }

  cat("min_score=", min_score, "\n");

  rp <- list( infp_bed = infp_bed, peak_broad=broadpeak_sum, min_score=min_score );

  rp <- start_calling( rp, min_score, pv_adjust, pv_threshold, smoothwidth, ncores )

  return(rp);
}

peak_calling<-function( asvm, gdm, bw_plus_path, bw_minus_path, infp_bed=NULL, use_rgtsvm=TRUE, min_score=NULL, pv_adjust="fdr", pv_threshold=0.05, smoothwidth=4, cpu_cores=1, gpu_cores=1 )
{
  if(!file.exists(bw_plus_path))
    stop( paste("Can't find the bigwig of plus strand(", bw_plus_path, ")"));

  if(!file.exists(bw_minus_path))
    stop( paste("Can't find the bigwig of minus strand(", bw_minus_path, ")"));

  #cat("[1]", as.character(Sys.time()), "\n");
  if(use_rgtsvm)
  {
    require(Rgtsvm);
    if( class(asvm)=="svm" && use_rgtsvm) class(asvm)<-"gtsvm";
    asvm <- Rgtsvm::predict.load( asvm, gpu_cores, verbose=T);
  }

  #cat("[2]", as.character(Sys.time()), "\n");
  if( is.null(infp_bed) )
  {
    infp_bed <- get_informative_positions(bw_plus_path, bw_minus_path, depth= 0, step=50, use_ANDOR=TRUE, use_OR=FALSE);
    infp_bed <- data.frame(infp_bed, pred=eval_reg_svm(gdm, asvm, infp_bed, bw_plus_path, bw_minus_path, batch_size= 50000, ncores=cpu_cores, use_rgtsvm=use_rgtsvm));
  }

  #cat("[3]", as.character(Sys.time()), "\n");
  colnames(infp_bed) <- c("chr", "start", "end", "pred");

  ## broad peaks with information are returned back.
  rp <- get_dense_infp( asvm, gdm, infp_bed, bw_plus_path, bw_minus_path, cpu_cores, use_rgtsvm);
  if( NROW(rp$infp_bed)>0 )
     colnames(rp$infp_bed) <- c("chr", "start", "end", "score", "infp");

  if(use_rgtsvm)
    Rgtsvm::predict.unload( asvm );

  #cat("[4]", as.character(Sys.time()), "\n");
  rp <- start_calling( rp, min_score, pv_adjust, pv_threshold, smoothwidth, cpu_cores )

  #cat("[5]", as.character(Sys.time()), "\n");
  return(rp);
}

start_calling<-function( rp, min_score, pv_adjust, pv_threshold, smoothwidth, ncores )
{
  if(!is.null(min_score))
    rp$min_score <- min_score;
  min_score <- rp$min_score;
  cat("min_score=", min_score, "\n");

  cor_mat <- build_cormat( rp$infp_bed, dist=20 );
  show(cor_mat);

  peak.idx <- which( rp$peak_broad$max>=min_score );
  rp$peak_broad <- rp$peak_broad[peak.idx,]; 

  #tmp.rdata = tempfile(".rdata");
  #save( rp, file=tmp.rdata);
  BLOCKWIDTH <- 1000;

  tmp.rdata.list <- c();
  for(chr in as.character(unique(rp$peak_broad$chr)) )
  {
     peak_broad <- rp$peak_broad[ as.character(rp$peak_broad$chr) == chr,]
     peak_broad <- peak_broad[order(peak_broad$start),,drop=F]
     k.sect <- 1:ceiling(NROW(peak_broad)/BLOCKWIDTH)

     for(k in k.sect)
     {
         rpx <- list( peak_broad=peak_broad, infp_bed=NULL, k=k);
         
         idx.min <- (k-1)*BLOCKWIDTH + 1;
         idx.max <- (k-1)*BLOCKWIDTH + BLOCKWIDTH;
         if(idx.max>NROW(peak_broad)) idx.max <- NROW(peak_broad);
         
         ki <- which( rp$infp_bed[,2] >= peak_broad[idx.min,]$start &
                rp$infp_bed[,3] <= peak_broad[idx.max,]$end )

         if( NROW(ki)>0 )
            rpx$infp_bed <- rp$infp_bed[ki,,drop=F]
       
         tmp.rdata = tempfile(".rdata");
         save(rpx, file=tmp.rdata);

         tmp.rdata.list <- c( tmp.rdata.list, tmp.rdata);
      }   
  }
  cat("block_number=", NROW(tmp.rdata.list), "\n");
 
  cpu.fun <- function( tmp.rdata )
  {
	require(dREG);
	
	#load 'rpx' object from splitted rdata file.
    load( tmp.rdata );
    if( is.null(rpx$infp_bed) ) return(NULL);
    unlink(tmp.rdata);
    
    file.RDS <- system.file("extdata", "rf-model-201803.RDS", package="dREG");
    model <- readRDS(file.RDS);

    idx.k <- (rpx$k-1)*BLOCKWIDTH + c(1:BLOCKWIDTH);
    idx.k <- idx.k[ idx.k <= NROW(rpx$peak_broad) ]

    P_list <- lapply(idx.k, function(kk){
      ki <- which( rpx$infp_bed[,2] >= rpx$peak_broad[kk,]$start &
                   rpx$infp_bed[,3] <= rpx$peak_broad[kk,]$end )
      xp <- rpx$infp_bed[ki,2]
      yp <- rpx$infp_bed[ki,4]
      if(NROW(xp) <= 3 || max(yp) <= min_score )
      {
        P <- NULL;
      }
      else
      {
		## Timestamp: 201803
        ## using Random Forest to split or merge the two adjacent local maximas.
        ## Removing any peaks less than 50
        P <- find_rf_peaks( model, xp, yp, SlopeThreshold=0.01, AmpThreshold=min_score, smoothwidth=smoothwidth, smoothtype=2, cor_mat=cor_mat);

		## Timestamp: 201712
        ## using 250 bp as the mimnial peak width to merge narrow peaks.
        ## Removing any peaks less than 100
        ## P <- find_peaks( xp, yp, SlopeThreshold=0.01, AmpThreshold=min_score, smoothwidth=smoothwidth, smoothtype=2, cor_mat=cor_mat);
      }

      if(!is.null(P)) P <- data.frame(chr=rpx$peak_broad[kk,]$chr, kk, P[,-1,drop=F]);
      return(P);
      });

    P_list <- as.data.frame(do.call("rbind", P_list))
    return(list(k=rpx$k, idx.k=idx.k, kk=idx.k, ret=P_list));
  }

  if(ncores>1)
  {
    sfInit(parallel = TRUE, cpus = ncores, type = "SOCK" )
    sfExport("BLOCKWIDTH", "min_score", "smoothwidth", "cor_mat" );

    fun <- as.function(cpu.fun);
    environment(fun)<-globalenv();

    dregP_list <- sfClusterApplyLB( tmp.rdata.list, fun=fun);
    sfStop();
  }
  else
    dregP_list <- lapply( tmp.rdata.list, cpu.fun );

  rp$raw_peak <- as.data.frame(do.call("rbind", lapply(dregP_list, function(x){if(!is.null(x)) return(x$ret) else return(NULL) })));
  colnames(rp$raw_peak) <- c("chr", "kk", "start", "end", "score", "prob", "smooth.mode", "original.mode", "centroid");

  if( NROW(rp$raw_peak)>0 )
    rp$peak_bed <- select_sig_peak( rp$raw_peak,  pv_adjust, pv_threshold );

  rp$peak_sum <- summary_peak( rp$raw_peak, rp$peak_bed );
  
  ## remove the "kk" column
  rp$raw_peak <- rp$raw_peak[, -2];

  return(rp);
}

select_sig_peak <- function( raw_peak, pv_adjust, pv_threshold )
{
  ## only select score, prob and center
  peak_bed <- raw_peak[,c("chr", "start", "end", "score", "prob", "original.mode"), drop=F];
  colnames(peak_bed) <- c("chr", "start", "end", "score", "prob", "center");

  ## remove peaks which width < 100bp;
  peak_bed <- peak_bed [ peak_bed$prob != -1,, drop=F];

  ## make multiple correction.
  peak_bed$prob <- p.adjust( peak_bed$prob, method = pv_adjust);

  ## only select significant peaks
  peak_bed <- peak_bed[ peak_bed$prob <= pv_threshold,, drop=F];

  return(peak_bed);
}


summary_peak <- function( raw_peak, peak_bed )
{
  peak_sum <- list();
  peak_sum$adjust.BH.0.05 <- sum( p.adjust(raw_peak$prob[ raw_peak$prob != -1 ], method="BH")<0.05)
  peak_sum$adjust.bonferroni.0.05 <- sum(p.adjust(raw_peak$prob[ raw_peak$prob != -1], method="bonferroni")<0.05)
  peak_sum$adjust.holm.0.05 <- sum(p.adjust(raw_peak$prob[ raw_peak$prob != -1], method="holm")<0.05)
  peak_sum$adjust.hochberg.0.05 <- sum(p.adjust(raw_peak$prob[ raw_peak$prob != -1], method="hochberg")<0.05)
  peak_sum$adjust.BY.0.05 <- sum(p.adjust(raw_peak$prob[ raw_peak$prob != -1], method="BY")<0.05)
  peak_sum$adjust.fdr.0.05 <- sum(p.adjust(raw_peak$prob[ raw_peak$prob != -1], method="fdr")<0.05)
  peak_sum$adjust.none.0.05 <- sum(p.adjust(raw_peak$prob[ raw_peak$prob != -1], method="none")<0.05)

  peak_sum$peak.narrow50 <- sum(raw_peak$prob == -1)

  peak <- cbind(raw_peak, prob2=raw_peak$prob)
  peak[which(peak$prob2!=-1),7]  <- p.adjust(peak$prob2[peak$prob2!=-1], method="fdr")
  peak_sum$peak.sig.score <- range(peak$score[which(peak[,7]<=0.05)])
  peak_sum$peak.narrow50.sig <- sum(raw_peak$prob==-1 & raw_peak$score > peak_sum$peak.sig.score [1] )
  peak_sum$peak.narrow50.score <- range(peak$score[peak$prob==-1] )

  return(peak_sum);
}


get_dense_infp <- function( asvm, gdm, infp_bed, bw_plus_path, bw_minus_path, ncores=1, use_rgtsvm=TRUE)
{
  pred_dense_infp<-function( dreg_peak, newinfp )
  {
	if(NROW(dreg_peak)==0) return(NULL);

    for(chr in unique( dreg_peak$chr ))
    {
      idx.chr <- which( as.character(dreg_peak$chr)==as.character(chr))
      if(NROW(idx.chr)==0) next;

      cat("=====", chr, "\n");

      predx <- dreg_peak[idx.chr, ];
      r.dense <- as.data.frame(rbindlist( lapply(1:NROW(predx), function(k){
        r.pos <- unique(c(seq(predx[k,]$start, predx[k,]$end, 10),  predx[k,]$end));
        return( data.frame( chr=chr, start=r.pos, end=r.pos+1));
      })));

      newinfp.chr <- newinfp[ newinfp$chr == as.character(chr),,drop=FALSE ];
      dup.rm <- which(paste(r.dense[,1], r.dense[,2], sep=":") %in% paste(newinfp.chr[,1], newinfp.chr[,2], sep=":"))
      if( NROW(dup.rm)>0 ) r.dense <- r.dense[-dup.rm,]

	  infp_dense <- data.frame(r.dense, pred=eval_reg_svm( gdm, asvm, r.dense, bw_plus_path, bw_minus_path, ncores=ncores, use_rgtsvm=use_rgtsvm, debug= TRUE));
      infp_dense<- infp_dense[ infp_dense$pred > 0.05,,drop=F];

      if(NROW(infp_dense)>0)
      {
        colnames(infp_dense) <- colnames(newinfp)[1:4];
        newinfp <- rbind( newinfp, cbind(infp_dense, INFP=0) );
      }
    }
    newinfp <- newinfp[with(newinfp, order(chr, start)),];
    return(newinfp);
  }

  ## get the min score for peak selection
  y.minus <- infp_bed$pred[ infp_bed$pred<0 ];
  laplace_sigma <- get_laplace_sigma( c(y.minus, abs(y.minus), infp_bed$pred[infp_bed$pred==0]))
  min_score <- get_laplace_quantile( laplace_sigma, 0.001 );

  ## fill the gaps between the two adjacent informative sites.
  gap_bed <- find_gap_infp( infp_bed, min_score, ncores=ncores );
  
  if(NROW(gap_bed)>0)
  {
    gap_score <- eval_reg_svm( gdm, asvm, gap_bed, bw_plus_path, bw_minus_path, ncores=ncores, use_rgtsvm=use_rgtsvm)
    gap_bed <- data.frame(gap_bed, pred=gap_score);

    ## adding the INFP to new informative sites.
    newinfp <- rbind( cbind(gap_bed, INFP=0), cbind(infp_bed, INFP=1));
  }
  else
    newinfp <- cbind(infp_bed, INFP=1);
  
  newinfp <- newinfp[with( newinfp, order(chr, start)),];

  broadpeak_sum <- get_broadpeak_summary( newinfp[,-5], threshold=0.05 );
  if (NROW(broadpeak_sum) != sum( !is.na(broadpeak_sum$max)) )
  {
    ## remove unknown contig, e.g.g chr1_gl000192_random chr10, chr17_ctg5_hap1, chr6_mann_hap4,chr6_qbl_hap6
    ## Notice: chr1_gl000192_random causes Bedmap is failed to get summary information for chr10-chr19
    ##         remove these contigs temporally
    infp_bed <- newinfp[ grep("_", newinfp$chr,invert=TRUE), ]
    broadpeak_sum <- get_broadpeak_summary( infp_bed[,-5], threshold=0.05 );
  }

  dense_infp <- pred_dense_infp( broadpeak_sum[ broadpeak_sum$max>=min_score,,drop=F ], newinfp );

  return(list(peak_broad=broadpeak_sum, infp_bed=dense_infp, min_score=min_score ));
}


find_gap_infp <- function( dreg_pred, threshold=0.2, ncores=1 )
{
  dreg_pred <- dreg_pred[with( dreg_pred, order(chr, start)),];

  cpu.fun<-function(chr)
  {
	require(data.table);

    predx <- dreg_pred[ as.character(dreg_pred$chr)==as.character(chr),];

    ## the distance between two ajacent informative sites.
    dist  <- predx[-1,]$start - predx[-NROW(predx),]$start;

    ## select the gap if distance > 50
    r.chr <- rbindlist( lapply(which(dist>50), function(k){
      r.pos <- c();

      ## if the current site has a high score
      if( predx[k,4] > threshold)
      {
        if(dist[k]<500)
          r.maxpos <- floor((dist[k]-50)/50)
        else
          r.maxpos <- ceiling((predx[k,4] - threshold)/0.05)

        r.pos <- c(r.pos, predx[k,2]+c(1:r.maxpos)*50);
      }

      ## if the neighbor site has a high score
      if( predx[k+1,4]>threshold)
      {
        if(dist[k]<500)
          r.maxpos <- floor((dist[k]-50)/50)
        else
          r.maxpos <- ceiling((predx[k+1,4] - threshold)/0.05)

        r.pos <- c(r.pos, predx[k+1,2]-c(r.maxpos:1)*50);
      }

      #remove duplicated postion
      r.pos <- unique(r.pos);

      if(NROW(r.pos)>0)
        return( data.frame( chr=chr, start=r.pos, end=r.pos+1))
      else
        return(c()); } ));

    return(r.chr);
  }

  if(ncores>1)
  {
    sfInit(parallel = TRUE, cpus = ncores, type = "SOCK" )
    sfExport("dreg_pred", "threshold" );

    fun <- as.function(cpu.fun);
    environment(fun)<-globalenv();

    r.gap <- sfClusterApplyLB( unique(dreg_pred$chr), fun=fun);
    sfStop();
  }
  else
    r.gap <- lapply( unique(dreg_pred$chr), cpu.fun );

  ## return a bed data containing all sites which are not predicted by the first run(informative sites)
  gap.bed <- as.data.frame(unique(rbindlist(r.gap)))
  if(NROW(gap.bed)==0)
     return(c());

  ## remove duplciated position
  dup.rm <- which(paste(gap.bed[,1], gap.bed[,2], sep=":") %in% paste(dreg_pred[,1], dreg_pred[,2], sep=":"))
  if( NROW(dup.rm)>0 ) gap.bed <- gap.bed[-dup.rm,]

  return( gap.bed);
}

get_broadpeak_summary <- function( infp_bed, threshold=0 )
{
  tb.peak <- merge_broad_peak(infp_bed, threshold);

  options("scipen"=100, "digits"=4)
  file.peak <- tempfile("peak", ".", ".peak.bed");
  write.table(data.frame( tb.peak[,-4 ], 1:NROW(tb.peak)), file.peak, col.names=F, row.names=F, quote=F, sep="\t");

  file.infp_bed <- tempfile("temp", ".", ".pred.bed");
  write.table(data.frame( infp_bed[,c(1,2,3)], "n", infp_bed[,4]), file.infp_bed, col.names=F, row.names=F, quote=F, sep="\t");

  tb.peak.sum <- read.table(pipe(paste("bedmap --echo --min --max --mean --sum --stdev --count --delim '\t' ", file.peak, file.infp_bed, sep=" ")), header=F);

  colnames(tb.peak.sum) <- c("chr", "start", "end", "no", "min", "max", "mean", "sum", "stdev","count");
  show(head(tb.peak.sum));

  unlink(c(file.peak, file.infp_bed));

  return(tb.peak.sum);
}


merge_broad_peak<-function( pred.bed, threshold, join = 500)
{
  ## bash command
  ##zcat $i | awk 'BEGIN{OFS="\t"} ($4 > '"$TH"') {print $1,$2-50,$3+51}' | sort-bed - | bedops --merge - | bed_merge.pl 500 > $j.no_score.bed

  pred.bed <- pred.bed[pred.bed[,4]>=threshold,,drop=F];
  if(NROW(pred.bed)==0) return(NULL);

  pred.bed[,2] <- pred.bed[,2] - 50;
  pred.bed[,3] <- pred.bed[,3] + 50;
  pred.bed <- pred.bed[with(pred.bed, order(chr, start)),];

  peak.bed <- c();
  for( chr in unique(as.character(pred.bed[,1])) )
  {
    pred.chr <- pred.bed[ as.character(pred.bed[,1]) == chr, ];

    dist <- pred.chr[ -1, 2 ] - pred.chr[ -NROW(pred.chr), 3 ];
    ## only select the long distance with right neighbor
    separe.idx <- which(dist >= 500);
    ## adding the first point
    separe.idx <- c(0, separe.idx);

    peak.idx <- data.frame(start=separe.idx[-NROW(separe.idx)]+1, end=separe.idx[-1])
    peak.idx <- peak.idx[ peak.idx[,1] != peak.idx[,2], ];

    if (NROW(peak.idx)>0)
    	peak.bed <- rbind(peak.bed,  data.frame(chr=chr, start=pred.chr[ peak.idx[,1] , 2 ], end=pred.chr[ peak.idx[,2], 3 ]) );
  }

  return( peak.bed );
}


find_peaks <- function( x, y, SlopeThreshold, AmpThreshold, smoothwidth, smoothtype=2, cor_mat=diag(5) )
{
  if (length(AmpThreshold)==1) AmpThreshold = rep( AmpThreshold, NROW(SlopeThreshold))
  if (length(smoothwidth)==1) smoothwidth = rep( smoothwidth, NROW(SlopeThreshold))
  #if (length(peakgroup)==1) peakgroup = rep( peakgroup, NROW(SlopeThreshold))

  smoothwidth <- round(smoothwidth);

  y.org <- y;
  y <- SegmentedSmooth( y, ifelse(NROW(y)>smoothwidth, smoothwidth, 3), smoothtype )

  if (smoothwidth>1)
    d <- SegmentedSmooth( deriv(y), ifelse(NROW(y)>smoothwidth, smoothwidth, 3), smoothtype )
  else
    d <- deriv(y);

  j <- ( 2 * round( smoothwidth/2)-1) : ifelse( (NROW(y) - smoothwidth-1)>0, (NROW(y) - smoothwidth-1),1);
  d.diff <- d[-1] * d[-NROW(d)];
  peak.loci <- which( d.diff <= 0 & y[-1]>AmpThreshold ) + 1;

  if(NROW(peak.loci)>1)
  {
	fake.rm.loci <- c();
    while( NROW(peak.loci)> 1 && any( x[peak.loci[-1]] - x[peak.loci[-NROW(peak.loci)]] <= 250) )
    {
      dist <- x[ peak.loci[-1] ] - x[ peak.loci[-NROW(peak.loci)]];

	  d.min <- which.min(dist);
	  if ( y[ peak.loci[ d.min ] ] >= y[ peak.loci[ d.min+1 ] ] && min(y[ peak.loci[d.min]:peak.loci[d.min+1] ]) < 0.5*y[ peak.loci[ d.min+1 ] ] )
	  {
		  fake.rm.loci <- c(fake.rm.loci, peak.loci [ d.min+1 ] );
          peak.loci <- peak.loci [ - (d.min+1) ]
      }
	  else if ( y[ peak.loci[ d.min ] ] < y [ peak.loci[ d.min+1 ] ] && min(y[ peak.loci[d.min]:peak.loci[d.min+1] ]) < 0.5*y[ peak.loci[ d.min ] ] )
	  {
		  fake.rm.loci <- c(fake.rm.loci, peak.loci [ d.min] );
          peak.loci <- peak.loci [ - d.min ]
      }
	  else if ( y[ peak.loci[ d.min ] ] < y[ peak.loci[ d.min+1 ] ] )
      {
        peak.rm <- peak.loci [ d.min ]
        peak.kp <- peak.loci [ d.min + 1  ]
        peak.loci <- peak.loci [ - d.min ]
        y[peak.rm:peak.kp] <- (y[peak.kp] - y[peak.rm] ) / ( peak.kp-peak.rm )*(c( peak.rm:peak.kp)-peak.rm) + y[peak.rm];
      }
      else
      {
        peak.rm <- peak.loci [ d.min +1]
        peak.kp <- peak.loci [ d.min ]
        peak.loci <- peak.loci [ -( d.min+1) ]
        y[peak.kp:peak.rm] <- y[peak.kp] - (y[peak.kp] - y[peak.rm] )/( peak.rm-peak.kp )*(c( peak.kp:peak.rm)-peak.kp) ;
      }

      next;
    }

    peak.loci <- sort( c(peak.loci, fake.rm.loci) );
  }

  P <- c(); i<-1;
  for( peak in peak.loci)
  {
    # the peak is at the left
    if(all(peak <= peak.loci))
    {
      i.left <- max(which( y[1:peak] <= AmpThreshold ))
      if(is.infinite(i.left))  i.left <- 1;
    }

    # the peak is at the right
    if(all(peak >= peak.loci))
    {
      i.right <- min(which( y[peak:NROW(y)] <= AmpThreshold )) + (peak-1)
      if(is.infinite(i.right))  i.right <- NROW(y);
    }

    if( any(peak.loci>peak) )
    {
      i.next <- peak.loci[min(which(peak.loci>peak))];
      suppressWarnings( i.valley0 <- min(which( y[peak:i.next]<AmpThreshold)) + peak -1)
      i.valley1 <- c(peak:i.next)[ which.min( y[peak:i.next] )];
      if(is.infinite(i.valley0))
        i.right <- i.valley1;
      if(!is.infinite(i.valley0) && !is.infinite(i.valley1))
      {
        if(i.valley0<i.valley1)
          i.right <- i.valley0
        else
          i.right <- i.valley1;
      }
    }

    if( any(peak.loci<peak) )
    {
      i.prev <- peak.loci[max(which(peak.loci<peak))];
      suppressWarnings( i.valley0 <- max(which( y[i.prev:peak]<AmpThreshold)) + i.prev -1)
      i.valley1 <- c(i.prev:peak)[which.min( y[i.prev:peak] ) ];
      if(is.infinite(i.valley0))
        i.left <- i.valley1 ;
      if(!is.infinite(i.valley0) && !is.infinite(i.valley1))
      {
        if(i.valley0>i.valley1)
          i.left <- i.valley0
        else
          i.left <- i.valley1;
      }
    }
    if( x[i.right] - x[i.left] >= 100 )
    {
      y.max <- max(y[i.left:i.right]);
      y.p <- which.max(y[i.left:i.right]) + i.left - 1;
      suppressWarnings( y.left <- min(which(y[i.left:y.p] >= y.max/5)) + i.left - 1 );
      if(is.infinite(y.left))  y.left<- i.left;
      suppressWarnings( y.right <- min(which(y[y.p:i.right] <= y.max/5)) + y.p - 1 );
      if(is.infinite(y.right)) y.right <- i.right;
      w.left <- x[y.p] - x[y.left];
      w.right <- x[y.right] - x[y.p];

      if(w.left > 2*w.right && w.right > 300 ) w.left <- 2*w.right;
      if(2*w.left < w.right && w.left  > 300 ) w.right <- 2*w.left;

      # Five Fields:  index, Peak, Height, left, right
      if( (x[y.p]-w.left+10) < ( x[y.p]+w.right-10) && (w.right+w.left>=100))
      {
        z.max <- -Inf;
        #the distace is 20 bp
        z.mi  <- NULL;
        for(i in 0:4)
        {
          z.i <- seq(y.p - i*2, y.p - i*2 + 8, 2);
          if(all(z.i>0) && all(z.i<=NROW(y.org)) && sum(y[z.i], na.rm=T)>z.max)
          {
            z.mi <- z.i
            z.max <- sum(y[z.i], na.rm=T);
          }
        }

        if( !is.null(z.mi) && sum(is.na(y[z.mi]))==0 )
        {
          pv <- pmvLaplace(y[z.mi], cor_mat);

#cat(pv, "==", y.org[z.mi], "\n");
          y.centriod <- sum( y[i.left:i.right]*c(1:(i.right-i.left+1)))/sum(y[i.left:i.right]);
          x.wc <- round(y.centriod/(i.right-i.left)*(x[i.right]-x[i.left])) + x[i.left];
          P <- rbind(P, c(i, x[y.p]-w.left+10, x[y.p]+w.right-10, max(y.org[i.left:i.right]), 1-pv, x[y.p], x[which.max(y.org[i.left:i.right])+i.left-1], x.wc) );
          i <- i+1;
        }
      }
    }
    else
    {
	   y.p <- which.max(y[i.left:i.right]) + i.left - 1;
       P <- rbind(P, c(i, x[i.left], x[i.right], max(y.org[i.left:i.right]), -1, x[y.p], 0, 0) );
    }

  }

  return(P);
}

build_cormat<-function(dreg_bed, dist=20, order=5)
{
  var.fun<-function(ypred, trunc=FALSE, cutoff=0.5)
  {
    y.suppose <- 0;
    #std <- sqrt( 2 * sum( abs(ypred - y.suppose)/NROW(ypred) ) ^2 );
    if(trunc)
    {
      outlier <-  which(abs(ypred - y.suppose) > cutoff );
      if(NROW(outlier)>0) ypred <- ypred[-outlier];
    }

    return(var(ypred - y.suppose, na.rm=T));
  }

  dist <- dreg_bed[-c(1,2),2] - dreg_bed[-(NROW(dreg_bed)-c(0,1)),2];
  cor.bed <- dreg_bed[ which(dist==20)+2,];
  cor.bed <- cor.bed[seq(1,NROW(cor.bed),2),]

  rho <- c(1)
  for(i in 1:4)
  {
    rho <- c(rho, cor(cor.bed$score[-c(1:i)], cor.bed$score[-NROW(cor.bed$score)+c(0:(i-1)) ]))
  }

  sigma2 <- var.fun( dreg_bed$score, trunc=T );
  mat <- matrix(1:order, nrow=order, ncol = order);
  cormat <-  sigma2*matrix( rho[abs( mat - t(mat) ) + 1 ], nrow=order, ncol=order);
  cormat;
}

get_laplace_sigma <- function( ypred, y=rep(0, NROW(ypred)) )
{
  na.idx <- which( is.na(ypred) | is.na(y) );
  if(NROW(na.idx)>0)
  {
    ypred <- ypred [ -na.idx];
    y     <- y [ -na.idx];
  }

  std <- sqrt( 2 * sum( abs(ypred - y)/NROW(y) ) ^2 );
  outlier <-  which(abs(ypred - y) > 5*std);
  if (length(outlier)>0)
    mae <- sum( abs(ypred - y)[-outlier])
  else
    mae <- sum( abs(ypred - y));

  mae  <- mae/( NROW(y)-length(outlier) );
  #cat ("Prob. model for test data: target value = predicted value + z,\nz: Laplace distribution e^(-|z|/sigma)/(2sigma),sigma=", mae, "\n" );
  return(mae);
}

get_laplace_quantile<-function(sigma, p=0.05)
{
  q <- qlaplace(1-p/2, m=0, s=sigma);
  return(q);
}

pmvLaplace<-function(x, cor_mat )
{
  ## when prob is approching 1, the different between normal distr. and laplace distr. is very few.
  p.norm <- pmvnorm(lower=-abs(x), upper=abs(x), mean=rep(0,length(x)), sigma=cor_mat)[1];
  if(p.norm>0.99)
    #return(c(p.norm, p.norm));
    return( p.norm );

  z <- c(10^-100, 10^seq(-20,-2, 1), seq( 0.02, 1, 0.04 ), seq(1, 10, 0.2), seq( 10.5, 100, 0.5 ), 10^(3:20), 10^100)
  p0 <- unlist( lapply(z, function(z0) {pmvnorm(lower=-abs(x/sqrt(z0)), upper=abs(x/sqrt(z0)), mean=rep(0,length(x)), sigma=cor_mat)*  exp(-z0)[1]} ) );
  p.max <- sum((z[-1] - z[-NROW(z)])*p0[-NROW(p0)]);
  if(p.max>1) p.max<-1
  p.min <- sum((z[-1] - z[-NROW(z)])*p0[-1]);
  if(p.min>1) p.min<-1

  #return(c(mean(p.max, p.min), p.norm) );
  return( mean(p.max, p.min) );
}



