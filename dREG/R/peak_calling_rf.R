find_rf_peaks <- function( model, x, y, SlopeThreshold, AmpThreshold, smoothwidth, smoothtype=2, cor_mat=diag(5) )
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
  #peak.loci <- which( d.diff < 0 & y[-1]>AmpThreshold ) + 1;
  peak.loci <- which( d[-1]<0 & d[-NROW(d)]> 0 & y[-1]>AmpThreshold ) + 1;

  if(NROW(peak.loci)==0)
     return(NULL);

  rp <- NULL;
  if(NROW(peak.loci)>1)
  {
	 dist <- x[ peak.loci[-1] ] - x[ peak.loci[-NROW(peak.loci)]];
	 LI <- peak.loci[-NROW( peak.loci)];
	 RI <- peak.loci[-1];
	 VI <- unlist( lapply(1:NROW(LI), function(i) c(LI[i]:RI[i])[which.min(y[LI[i]:RI[i]])] ) );

	 start <- x[LI];
	 stop <- x[RI];
	 valley <- x[VI]

	 LS <- y[LI];
	 RS <- y[RI];
	 VS <- y[VI];

     rp <- data.frame(dist=dist, LI=LI, RI=RI, VI=VI, start=start, stop=stop, valley=valley, LS=LS, RS=RS, VS=VS, ST=0);
  }

  #the most left side
  rp.left <- c(dist=x[peak.loci[1]]-x[1], LI=1, RI=peak.loci[1], VI=1,
                   start=x[1], stop=x[peak.loci[1]], valley =x[1],
                   LS=y[1], RS=y[peak.loci[1]], VS=y[1], ST=-1);
  #the most right side
  rp.right <- c(dist=x[NROW(x)]-x[peak.loci[NROW(peak.loci)]], LI=peak.loci[NROW(peak.loci)], RI=NROW(x), VI=NROW(x),
                   start=x[peak.loci[NROW(peak.loci)]], stop=x[NROW(x)], valley=x[NROW(x)],
                   LS=y[peak.loci[NROW(peak.loci)]], RS=y[NROW(x)], VS=y[NROW(x)], ST=1);

  rp <- rbind(rp.left, rp, rp.right);

  rp0 <- rp;
  rp <- split_peak( model, rp );
  if(NROW(rp)==0)
  {
    show(rp0);
    browser();
    return(NULL);
  }

  P <- c();
  for( i in 1:NROW(rp) )
  {
    i.left <-  min(which ( y[rp[i,"LI"]:rp[i,"RI"]] > AmpThreshold) );
    i.right <- max(which ( y[rp[i,"LI"]:rp[i,"RI"]] > AmpThreshold) );
    i.peak <-  which.max ( y[rp[i,"LI"]:rp[i,"RI"]]);
    if( NROW(i.left)==0 || NROW(i.right)==0 || NROW(i.peak)==0 )
       next;

    i.left <-  c(rp[i,"LI"]:rp[i,"RI"])[i.left];
    i.right <-  c(rp[i,"LI"]:rp[i,"RI"])[i.right];
    i.peak <-  c(rp[i,"LI"]:rp[i,"RI"])[i.peak];

	w.left <- x[i.peak] - x[i.left];
	w.right <- x[i.right] - x[i.peak];

    if(w.left > 2*w.right && w.right > 300 ) w.left <- 2*w.right;
    if(2*w.left < w.right && w.left  > 300 ) w.right <- 2*w.left;

    if( i.right - i.left < 5 )
    {
	  y.p <- which.max(y[i.left:i.right]) + i.left - 1;
	  P <- rbind(P, c(i, x[i.left], x[i.right], max(y.org[i.left:i.right]), -1, x[y.p], 0, 0) );
    }
    else
	{
	  ## the failure case: 3<->11 ,if peak=6, i.sample=4,6,8,10,12
      if( i.right - i.left < 9)
        i.sample <- sort(c(i.left:i.right)[order(y[i.left:i.right], decreasing=T)][1:5])
      else
      {
         mat.sample <- rbind(c(i.peak-4, i.peak-2, i.peak-0, i.peak+2, i.peak+4),
                             c(i.peak-2, i.peak-0, i.peak+2, i.peak+4, i.peak+6),
                             c(i.peak-6, i.peak-4, i.peak-2, i.peak-0, i.peak+2),
                             c(i.peak-8, i.peak-6, i.peak-4, i.peak-2, i.peak+0),
                             c(i.peak-0, i.peak+2, i.peak+4, i.peak+6, i.peak+8));

        i.true <- which( unlist(apply(mat.sample, 1, function(x) { all(x %in% c(i.left:i.right) ) } )) ) [1];
        i.sample <- mat.sample[i.true,]

        if(is.na(i.true))
        {
           show(mat.sample);
           cat(i.left, i.right, i.peak, i.true, "\n");
           cat(i.sample, "\n");
        }
      }

cat(y[i.sample], "\n\n");

      pv <- NA;
      if( !is.null(i.sample) && sum(is.na(y[i.sample]))==0 )
      pv <- pmvLaplace(y[i.sample], cor_mat);

      y.centriod <- sum( y[i.left:i.right]*c(1:(i.right-i.left+1)))/sum(y[i.left:i.right]);
      x.wc <- round(y.centriod/(i.right-i.left)*(x[i.right]-x[i.left])) + x[i.left];
      P <- rbind(P, c( i, x[i.peak]-w.left+10, x[i.peak]+w.right-10, max(y.org[i.left:i.right]), 1-pv, x[i.peak],
                 x[which.max(y.org[i.left:i.right])+i.left-1], x.wc) );
    }
  }

  return(P);
}


##
## rp format:
##
##
##start :
##stop  :
##vally :
##   LI : left index
##   RI : right index
##   VI : valley index
##   LS : left maximum score
##   RS : right maximum  score
##   VS : peak valley score
##   ST : -2: removed, -1(left side), 0(unknown), 1(right side), 2:merge, 3: split
##  dist: distance(LD+RD)
##(*)LD : Left distance
##(*)RD : right distance
##(*)MS : max score, max(LS,RS)
##(*)d1 : max(LS,RS)-min(LS,RS)
##(*)d2 : min(LS,RS)-vs
##(*)dr : d2/(d1+VS)

split_peak<-function(model, rp )
{
  rowMins <-function(x){ apply(x, 1, min) }
  rowMaxs <-function(x){ apply(x, 1, max) }

  rp <- data.frame(IDX=1:NROW(rp), rp);
  rp$LD <- rp$valley-rp$start;
  rp$RD <- rp$stop-rp$valley;
  rp$maxy <- rowMaxs(rp[,c("LS", "RS")]);
  rp$d1 <- abs(rp$LS-rp$RS)
  rp$d2 <- rowMins(rp[,c("LS", "RS")]) - rp$VS;
  rp$d3 <- rp$VS;
  rp$dr <- rp$d2/(rp$d1+rp$VS);

  while(sum(rp$ST==0)>0)
  {
    newdata <- rp[rp$ST==0, c("dist","LD", "RD", "LS","RS","maxy","d1","d2","d3", "dr")]
    colnames(newdata) <- c("dist", "r1", "r2", "y1", "y2", "maxy", "d1", "d2", "d3", "dr")
    pred <- predict(model, newdata=newdata);
    rp$ST[rp$ST==0] = c(2,3)[as.numeric(pred>0.5)+1];

    # search successive merge(3) region to do merge
    rcomp <- cbind( rp[-NROW(rp), c("IDX", "ST")], rp[-1, c("IDX","ST")] );
    colnames(rcomp) <- c("IDX", "ST0", "IDX1", "ST1");
    rcomp <- rcomp[rcomp$ST0==2 & rcomp$ST1==2 & rcomp$IDX+1==rcomp$IDX1, ];
    if(NROW(rcomp)==0) break;

    # only get the top idx by each region.
    idx   <- which(rp$IDX==rcomp$IDX[1]);
    idx1   <- idx + 1
    rp[idx+1,"ST"] <- -2;
    rp[idx,  "ST"] <- 0;
    rp[idx,  "stop"] <- rp[idx1,  "stop"];
    rp[idx,  "RI"] <- rp[idx1,"RI"];
    rp[idx,  "RS"] <- rp[idx1,"RS"];
    rp[idx,  "VI"] <- c(rp[idx1,"VI"], rp[idx,"VI"])[which.min( c(rp[idx1,"VS"], rp[idx,"VS"])) ];
    rp[idx,  "valley"]  <- c(rp[idx1,"valley"], rp[idx,"valley"])[which.min( c(rp[idx1,"VS"], rp[idx,"VS"])) ];
    rp[idx,  "LD"] <- rp[idx,  "valley"] - rp[idx,  "start"];
    rp[idx,  "RD"] <- rp[idx,  "stop"] - rp[idx,  "valley"]  ;
    rp[idx,  "dist"] <- rp[idx,  "stop"] - rp[idx,  "start"];
    rp[idx,  "VS"] <- min(rp[idx1,"VS"], rp[idx,"VS"]);
    rp[idx,  "d1"] <- max(rp[idx,"LS"], rp[idx,"RS"]) - min(rp[idx,"LS"], rp[idx,"RS"]);
    rp[idx,  "d2"] <- min(rp[idx,"LS"], rp[idx,"RS"]) - rp[idx,  "VS"];
    rp[idx,  "d3"] <- rp[idx,"VS"];
    rp[idx,  "dr"] <- rp[idx, "d2"]/( rp[idx,"d1"]+rp[idx, "d3"] );
    rp[idx,  "maxy"] <- max(rp[idx, "LS"], rp[idx, "RS"]);

    ## remove all merged regions
    rp <- rp[rp$ST != -2,];
    rp$IDX <- 1:NROW(rp);
  }

  # split regions
  rp$IDX <- rp$IDX*2;
  idx <- which( rp$ST==3 );
  if(NROW(idx)>0)
  {
    newr <- rbindlist( lapply(idx, function(i){
          rRegion <- rbind(rp[i, ], rp[i, ]);
          rRegion[1, "ST"] <- +1;
          rRegion[1, "stop"] <- rp[i,  "valley"];
          rRegion[1, "RS"] <- rp[i,  "VS"];
          rRegion[1, "RI"]  <- rp[i,"VI"]

          rRegion[2, "ST"] <- -1;
          rRegion[2, "IDX"]  <- rp[i,"IDX"]+1
          rRegion[2, "LI"]  <- rp[i,"VI"]
          rRegion[2, "LS"] <- rp[i,  "VS"];
          rRegion[2, "start"] <- rp[i, "valley"];
          return(rRegion)}));
    rp <- rbind( rp[-idx,], newr);
    rp <- rp[order(rp$IDX),];
  }

  rp <- as.data.frame(rp);
  rpeak <- data.frame(matrix(ncol = 9, nrow = 0))
  #make peak based on status=-1 and status=1
  LI <- RI <- PI <- start <- stop <- peak <- LS <- RS <- PS <- center <- center.s <- c();
  for(i in 1:NROW(rp))
  {
    if(rp[i, "ST"]==-1)
    {
      LI    <- rp[i, "LI"];
      LS    <- rp[i, "LS"];
      start <- rp[i, "start"];
      center <- c(rp[i, "start"], rp[i, "stop"]);
      center.s <- c(rp[i, "LS"], rp[i, "RS"]);
      center.i <- c(rp[i, "LI"], rp[i, "RI"]);
    }
    else if (rp[i, "ST"]==2)
    {
      center <- c(center, rp[i, "start"], rp[i, "stop"], rp[i, "vally"]);
      center.s <- c(center.s, rp[i, "LS"], rp[i, "RS"],rp[i, "VS"]);
      center.i <- c(center.i, rp[i, "LI"], rp[i, "RI"],rp[i, "VI"]);
    }
    else if (rp[i, "ST"]==1)
    {
      RI   <- rp[i, "RI"];
      RS   <- rp[i, "RS"];
      stop <- rp[i, "stop"];
      center <- c(center, rp[i, "start"], rp[i, "stop"]);
      center.s <- c(center.s, rp[i, "LS"], rp[i, "RS"]);
      center.i <- c(center.i, rp[i, "LI"], rp[i, "RI"]);

      peak <- center[which.max(center.s)];
      PI<- center.i[which.max(center.s)];
      rpeak <- rbind(rpeak, c(LI=LI, RI=RI, PI=PI, start=start, stop=stop, center=peak, LS=LS, RS=RS, VS=max(center.s)));
      LI <- RI <- PI <- start <- stop <- peak <- LS <- RS <- VS <- center <- center.s <- center.i <- NULL;
    }
  }

  colnames(rpeak) <- c("LI", "RI", "PI", "start", "stop", "center", "LS", "RS", "PS");
  return(rpeak);
}