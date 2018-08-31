check_bigwig<-function(bw_path, strand="+", out.file="")
{
  if(!file.exists(bw_path))
	stop( paste("Can't find the bigwig of plus strand(", bw_path, ")"));

  ## Load bigWigs
  bw  <- load.bigWig(bw_path)
  b.Norm <- b.Strand <- b.pointReads <- TRUE
  per.pointReads <- c();

  q_chroms <- bw$chroms[bw$chromSizes > 100*1000]
  for(chr in q_chroms)
  {
    chr_n <- which( bw$chroms == chr )
    pos <- round(runif(1, 1, bw$chromSizes[chr_n]-100*1000));

    qbw <- step.bpQuery.bigWig(bw, bw$chroms[chr_n], pos, min( pos + 1000000,bw$chromSizes[chr_n]), 1, with.attributes=FALSE);

    # Normalization checking
    b.Norm <- b.Norm & all(round(qbw)==qbw);

    ## not allow mixture of negative and positive values
  	b.Strand <- b.Strand & ( all(qbw>=0) || all(qbw<=0) );

    # Checking reads mapped at a point or a region
	idx <- which( qbw == 1 );
	if( length(idx)==0 )
		idx <- which( qbw == -1 );

    if( length(idx)>0 )
    {
		i.start <- i.stop <- idx[1];
		df.reg <- c();
    	for (i in idx)
    	{
			if (i.stop == i - 1 )
				i.stop <- i
			else
			{
				df.reg <- rbind(df.reg, c(i.start, i.stop));
				i.start <- i;
				i.stop  <- i;
			}
		}

		if( NROW(df.reg) > 0 )
		{
			distance <- df.reg[,2] - df.reg[,1]
			## at least 20% is the locus(point), not a region
			if( NROW(distance)>=100 )
				per.pointReads  <- c( per.pointReads , sum(distance==0)/NROW(distance) );
		}
	}
  }

  unload.bigWig(bw);

  if( NROW(per.pointReads)>1 )
	  b.pointReads <- (mean(per.pointReads)>0.8);

#cat("mean(per.pointReads", mean(per.pointReads), "\n");

  if(!b.Norm)
	 cat("The bigwig file might be normalized.\n", file=out.file, append=TRUE);

  if(!b.Strand)
	 cat("Read count should be >= 0 in plus strand and <= 0 in minus strand.\n", file=out.file, append=TRUE);

  if(!b.pointReads)
	 cat("Every read might be mapped to a region, not a locus.\n", file=out.file, append=TRUE);

  return( b.Norm & b.Strand & b.pointReads );
}