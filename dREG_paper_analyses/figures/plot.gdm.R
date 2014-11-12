## 
## Uses the grid package to plot the signal returned in a GDM.
plot.gdm <- function(gdm, positive_positions= NULL, bigwig_plus= NULL, bigwig_minus= NULL, gdata= NULL) {
  require(grid) ## To execute the visualization ... two more dependencies I'd rather not enforse on install.
  require(boot)

  stopifnot(!(is.null(positive_positions) & is.null(bigwig_plus) & is.null(bigwig_minus)) | !is.null(gdata))  ## At least one of these combos must be speicifed.
#  stopifnot(xor(!(is.null(positive_positions) & is.null(bigwig_plus) & is.null(bigwig_minus)), !is.null(gdata)))  ## Only specify one.

  ## Read logicstic representation of genomic data.
  if(is.null(gdata)) {
    gdata <- read_genomic_data(gdm, positive_positions, bigwig_plus, bigwig_minus, as_matrix= TRUE)
  }
  
  if(NROW(gdata) > 1) {
    dBoot <- boot(gdata, function(a, i) {colSums(a[i,])/ n_eval}, R=10000)
    std.error <- sapply(1:NCOL(gdata) , function(x) {sd(dBoot$t[,x])})
  }
  if(NROW(gdata) == 1) {
    dBoot <- list()
    dBoot$t0 <- gdata
	std.error <- rep(0, NROW(gdata))
  }
  
  fwd_sig <- c(1:(NCOL(gdata)/2))
  rev_sig <- c((NCOL(gdata)/2 +1):NCOL(gdata))

  drawBars <- function(signal, x_scale, half_nWindows, window_sizes, min.at.0=FALSE) {
    if(!min.at.0) {
      plus_labs= pretty(c(signal[["plus"]]), n=3)
      minus_labs= pretty(c(signal[["minus"]]), n=3)
	  range_plus <- max(plus_labs) - min(plus_labs)
	  range_minus <- max(minus_labs) - min(minus_labs)

      signal[["plus"]] <- signal[["plus"]] - min(plus_labs)
      signal[["minus"]] <- signal[["minus"]] - min(minus_labs)
    }
	else {
	  print("WARNING!  NO SUPPORT YET!")
      plus_labs= pretty(c(0, signal[["plus"]], signal[["minus"]]))
      minus_labs= pretty(c(0, signal[["plus"]], signal[["minus"]]))	
	}
  
    ## Print the barplot.
    next.vp <- viewport(x= 0.1, y= 0.5, width = 0.9, height = 1, just=c("left","center"))
    pushViewport(next.vp)

	trans<-half_nWindows*window_sizes 
    normx <- 2*trans
	width <- (window_sizes/3)/normx
	centers <- (c(seq(-half_nWindows*window_sizes,-window_sizes,window_sizes)+(window_sizes/4),
	           seq(window_sizes,half_nWindows*window_sizes, window_sizes)-(window_sizes/4)) +trans)/normx
	
    for(i in 1:length(centers)) {
      x <- c(centers[i]+width, centers[i]-width, centers[i]-width, centers[i]+width)
      yp <- c(0.55, 0.55, 0.55+c(signal[["plus"]][i], signal[["plus"]][i])/range_plus/2)
      ym <- c(0.45, 0.45, 0.45-c(signal[["minus"]][i], signal[["minus"]][i])/range_minus/2)
	 
 	  ## Draw box.
      grid.polygon(x, yp, gp=gpar(fill="#c5000b"))
      grid.polygon(x, ym, gp=gpar(fill="#0084d1"))
	 
      ## Draw std. error.
      grid.lines(x=rep(mean(x),2), y=0.55+c(signal[["plus"]][i], signal[["plus"]][i]+signal[["err.plus"]][i])/range_plus/2)
      grid.lines(x=x[1:2], y=0.55+rep(signal[["plus"]][i]+signal[["err.plus"]][i], 2)/range_plus/2)
      grid.lines(x=rep(mean(x),2), y=0.45-c(signal[["minus"]][i], signal[["minus"]][i]+signal[["err.minus"]][i])/range_minus/2)
      grid.lines(x=x[1:2], y=0.45-rep(signal[["minus"]][i]+signal[["err.minus"]][i], 2)/range_minus/2)
    }
	popViewport()
   
    ## Print a basic axis.
    next.vp <- viewport(x= 0.07, y= 0.5, width = 0.03, height = 1, just=c("left","center"))
	pushViewport(next.vp)
	  grid.yaxis(at= seq(0.55, 1, length= length(plus_labs)), label= plus_labs)
	  grid.yaxis(at= seq(0, 0.45, length= length(minus_labs)), label= rev(minus_labs))
	popViewport()
  }

  
  ## Prep the grid object.
  grid.newpage()
  top.vp <- viewport(width = 0.98, height = 0.98, xscale= c(0,1), yscale= c(0,1))
  pushViewport(top.vp)

  ## Get a matrix of the actual values.
  signal <- list()
  for(i in 1:gdm@n_zooms) {
    indx <- (2*sum(gdm@half_nWindows[0:(i-1)])+1):sum(gdm@half_nWindows[0:(i)]*2)
	signal[[i]] <- list()
    signal[[i]][["plus"]] <- dBoot$t0[fwd_sig][indx]
    signal[[i]][["minus"]] <- dBoot$t0[rev_sig][indx]
    signal[[i]][["err.plus"]] <- std.error[fwd_sig][indx]
    signal[[i]][["err.minus"]] <- std.error[rev_sig][indx]
	
    ## Plots bootstrapped data.
    next.vp <- viewport(width = 1, height = 0.9/gdm@n_zooms, y= (gdm@n_zooms-i)*1/gdm@n_zooms, xscale= c(0,1), yscale= c(0,1), just=c("center","bottom"))
    pushViewport(next.vp)
      drawBars(signal[[i]], half_nWindows=gdm@half_nWindows[i], window_sizes=gdm@window_sizes[i])
	popViewport()
  }

  popViewport()
}

