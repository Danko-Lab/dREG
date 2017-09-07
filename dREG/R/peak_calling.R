peak_calling<-function( asvm, gdm, bw_plus_path, bw_minus_path, infp_bed=NULL, ncores=1, use_rgtsvm=TRUE, min_score=NULL, smoothwidth=4 )
{
  if( is.null(infp_bed) )
  {
    infp_bed <- get_informative_positions(bw_plus_path, bw_minus_path, depth= 0, step=50, use_ANDOR=TRUE, use_OR=FALSE);
    infp_bed <- data.frame(infp_bed, pred=eval_reg_svm(gdm, asvm, infp_bed, bw_plus_path, bw_minus_path, batch_size= 50000, ncores=ncores, use_rgtsvm=use_rgtsvm));
  }

  colnames(infp_bed) <- c("chr", "start", "end", "pred");
  rp <- get_dense_infp( asvm, gdm, infp_bed, bw_plus_path, bw_minus_path, ncores, use_rgtsvm);

  cor_mat <- build_cormat( rp$infp_bed, dist=20 );
  show(cor_mat);

  if(is.null(min_score))
     min_score <- rp$min_score;

  cat("min_score=", min_score, "\n");	
  
  peak.idx <- which( rp$peak_sum$max>=min_score );
  k.sect <- 1:ceiling(NROW(peak.idx)/500)

  cpu.fun <- function(k)
  {
	require(dREG);
    load(tmp.rdata);

    idx.k <- (k-1)*500 + c(1:500);
    idx.k <- idx.k[ idx.k <= NROW(peak.idx) ]

    P_list <- lapply(peak.idx[idx.k], function(kk){
      ki <- which( as.character(rp$infp_bed[,1])== as.character(rp$peak_sum[kk,]$chr) &
                                rp$infp_bed[,2] >= rp$peak_sum[kk,]$start &
                                rp$infp_bed[,3] <= rp$peak_sum[kk,]$end )
      xp <- rp$infp_bed[ki,2]
      yp <- rp$infp_bed[ki,4]
      if(NROW(xp)<=3 || max(yp)<=0.1 )
      {
        P <- NULL;
      }
      else
      {
        #P <- try( find_peaks( xp, yp, SlopeThreshold=0.01, AmpThreshold=min_score, smoothwidth=smoothwidth, smoothtype=2, cor_mat=cor_mat) );
        P <- find_peaks( xp, yp, SlopeThreshold=0.01, AmpThreshold=min_score, smoothwidth=smoothwidth, smoothtype=2, cor_mat=cor_mat);
        if(class(P)=="try-error")
          P<-NULL;
      }

      if(!is.null(P)) P <- data.frame(chr=rp$peak_sum[kk,]$chr, kk, P[,-1,drop=F]);
      return(P);
        });

    P_list <- as.data.frame(do.call("rbind", P_list))
    return(list(k=k, idx.k=idx.k, kk=peak.idx[idx.k], ret=P_list));
  }

  tmp.rdata = tempfile(".rdata");
  save( rp, file=tmp.rdata);

  if(ncores>1)
  {
    sfInit(parallel = TRUE, cpus = ncores, type = "SOCK" )
    sfExport("tmp.rdata", "min_score", "smoothwidth", "cor_mat", "peak.idx");

    fun <- as.function(cpu.fun);
    environment(fun)<-globalenv();

    dregP_list <- sfClusterApplyLB( k.sect, fun=fun);
    sfStop();
  }
  else
    dregP_list <- lapply(  k.sect, cpu.fun );

  unlink(tmp.rdata);

  dregP <- as.data.frame(do.call("rbind", lapply(dregP_list, function(x){if(!is.null(x)) return(x$ret) else return(NULL) })));
  dregP <- dregP [,-2, drop=F];
  colnames(dregP) <- c("chr", "start", "end", "score", "prob.ml", "prob.mn", "smooth.mode", "original.mode", "centroid");

  rp$peak_bed <- dregP[dregP$prob_ml<=0.05,, drop=F];
  rp$peak_sum <- rp$peak_sum[rp$peak_sum$max>=min_score,,drop=F];

  return(rp);
}

get_dense_infp <- function( asvm, gdm, infp_bed, bw_plus_path, bw_minus_path, ncores=1, use_rgtsvm=TRUE)
{
  pred_dense_infp<-function( dreg_peak, newinfp )
  {
	if(NROW(dreg_peak)==0) return(NULL);

    for(chr in unique( dreg_peak$chr ))
    {
      cat("=====", chr, "\n");
      idx.chr <- which( as.character(dreg_peak$chr)==as.character(chr))
      if(NROW(idx.chr)==0) next;

      predx <- dreg_peak[idx.chr, ];
      r.dense <- as.data.frame(rbindlist( lapply(1:NROW(predx), function(k){
        r.pos <- unique(c(seq(predx[k,]$start, predx[k,]$end, 10),  predx[k,]$end));
        return( data.frame( chr=chr, start=r.pos, end=r.pos+1));
      })));

      dup.rm <- which(paste(r.dense[,1], r.dense[,2], sep=":") %in% paste(newinfp[,1], newinfp[,2], sep=":"))
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
  gap_score <- eval_reg_svm( gdm, asvm, gap_bed, bw_plus_path, bw_minus_path, ncores=ncores, use_rgtsvm=use_rgtsvm)
  gap_bed <- data.frame(gap_bed, pred=gap_score);

  ## adding the INFP to new informative sites.
  newinfp <- rbind( cbind(gap_bed, INFP=0), cbind(infp_bed, INFP=1));
  newinfp <- newinfp[with( newinfp, order(chr, start)),];

  peak_sum <- get_peak_summary( newinfp[,-5], threshold=0.05 );
  dense_infp <- pred_dense_infp( peak_sum[ peak_sum$max>=min_score,,drop=F ], newinfp );

  return(list(peak_sum=peak_sum, infp_bed=dense_infp, min_score=min_score ));
}


find_gap_infp <- function( dreg_pred, threshold=0.2, ncores=1 )
{
  dreg_pred <- dreg_pred[with( dreg_pred, order(chr, start)),];

  r.gap <- list();
  for(chr in unique(dreg_pred$chr))
  {
    predx <- dreg_pred[ as.character(dreg_pred$chr)==as.character(chr),];

    ## the distance between two ajacent informative sites.
    dist  <- predx[-1,]$start - predx[-NROW(predx),]$start;

    ## select the gap if distance > 50
    r.gap[[chr]] <- rbindlist( mclapply(which(dist>50), function(k){
      r.pos <- c();

      ## if thecurrent site has a high score
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

      if(NROW(r.pos)>0)
        return( data.frame( chr=chr, start=r.pos, end=r.pos+1))
      else
        return(c()); }, mc.cores = ncores ));
  }

  ## return a bed data containing all sites which are not predicted by the first run(informative sites)
  return( as.data.frame(unique(rbindlist(r.gap))));
}

get_peak_summary <- function( infp_bed, threshold=0 )
{
  tb.peak <- merge_broad_peak(infp_bed, threshold);

  options("scipen"=100, "digits"=4)
  file.peak <- tempfile("peak", ".", "bed");
  write.table(data.frame( tb.peak[,-4 ], 1:NROW(tb.peak)), file.peak, col.names=F, row.names=F, quote=F, sep="\t");

  file.infp_bed <- tempfile("temp", ".", "pred.bed");
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
  for( chr in unique(pred.bed[,1]) )
  {
    pred.chr <- pred.bed[ pred.bed[,1] == chr, ];

    dist <- pred.chr[ -1, 2 ] - pred.chr[ -NROW(pred.chr), 3 ];
    ## only select the long distance with right neighbor
    separe.idx <- which(dist >= 500);
    ## adding the first point
    separe.idx <- c(0, separe.idx);

    peak.idx <- data.frame(start=separe.idx[-NROW(separe.idx)]+1, end=separe.idx[-1])
    peak.idx <- peak.idx[ peak.idx[,1] != peak.idx[,2], ];
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
    while( NROW(peak.loci)> 1 && any( x[peak.loci[-1]] - x[peak.loci[-NROW(peak.loci)]] <= 250) )
    {
      dist <- x[ peak.loci[-1] ] - x[ peak.loci[-NROW(peak.loci)]];

      if ( y[ peak.loci[ which.min(dist) ] ] < y[ peak.loci[ which.min(dist)+1 ] ] )
      {
        peak.rm <- peak.loci [ which.min(dist) ]
        peak.kp <- peak.loci [ which.min(dist) + 1  ]
        peak.loci <- peak.loci [ - which.min(dist) ]
        y[peak.rm:peak.kp] <- (y[peak.kp] - y[peak.rm] ) / ( peak.kp-peak.rm )*(c( peak.rm:peak.kp)-peak.rm) + y[peak.rm];
      }
      else
      {
        peak.rm <- peak.loci [ which.min(dist) +1]
        peak.kp <- peak.loci [ which.min(dist) ]
        peak.loci <- peak.loci [ -( which.min(dist)+1) ]
        y[peak.kp:peak.rm] <- y[peak.kp] - (y[peak.kp] - y[peak.rm] )/( peak.rm-peak.kp )*(c( peak.kp:peak.rm)-peak.kp) ;
      }
      next;
    }
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
    rho <- c(rho, cor(cor.bed$pred[-c(1:i)], cor.bed$pred[-NROW(cor.bed$pred)+c(0:(i-1)) ]))
  }

  sigma2 <- var.fun( dreg_bed$pred, trunc=T );
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
    return(c(p.norm, p.norm));

  z <- c(10^-100, 10^seq(-20,-2, 1), seq( 0.02, 1, 0.04 ), seq(1, 10, 0.2), seq( 10.5, 100, 0.5 ), 10^(3:20), 10^100)
  p0 <- unlist( lapply(z, function(z0) {pmvnorm(lower=-abs(x/sqrt(z0)), upper=abs(x/sqrt(z0)), mean=rep(0,length(x)), sigma=cor_mat)*  exp(-z0)[1]} ) );
  p.max <- sum((z[-1] - z[-NROW(z)])*p0[-NROW(p0)]);
  if(p.max>1) p.max<-1
  p.min <- sum((z[-1] - z[-NROW(z)])*p0[-1]);
  if(p.min>1) p.min<-1

  return(c(mean(p.max, p.min), p.norm) );
}



