  drawBars <- function(signal, error, label) {
    ## Remove NAs.
    indx <- !is.na(signal)
    signal <- signal[indx]
    error <- error[indx]
    label <- label[indx]

    ## compute other parts...
    lab<- pretty(signal, n=5)
    yrange <- max(lab, na.rm=TRUE) - min(lab, na.rm=TRUE)
    width <- 1/NROW(signal)/3 
    centers <- seq(width, 1, 1/NROW(signal))
    zero <- 1-max(lab, na.rm=TRUE)/yrange

    ord <- order(signal)
    signal <- signal[ord]
    error  <- error[ord]
    label  <- label[ord]
	
    ## Print the barplot.
    require(grid)
    grid.newpage()
    top.vp <- viewport(width = 0.98, height = 0.98, xscale= c(0,1), yscale= c(0,1))
    pushViewport(top.vp)
    next.vp <- viewport(x= 0.1, y= 0.15, width = 0.9, height = 0.85, just=c("left","bottom"))
    pushViewport(next.vp)
    for(i in lab) {
      grid.lines(y= (i-min(lab))/yrange, gp=gpar(col="light gray"))
    }
  
    for(i in 1:NROW(signal)) {
      x <- c(centers[i]+width, centers[i]-width, centers[i]-width, centers[i]+width)
      ym <- signal[i]-min(lab)
      y <- c(zero, zero, c(ym, ym)/yrange)
	  
      if(zero>(ym/yrange)) {
         fill="light blue"
         ye <- ym-error[i]
      }
      else {
         fill="dark red"
         ye <- ym+error[i]
      }

      grid.polygon(x, y, gp=gpar(fill=fill))
	
   ## Draw std. error.
      grid.lines(x=centers[i], y=c(ym, ye)/yrange)
      grid.lines(x=x[1:2], y=(ye)/yrange)
    }
    popViewport()

    ## Plot Y axis.
    next.vp <- viewport(x= 0.07, y= 0.15, width = 0.03, height = 0.85, just=c("left","bottom"))
    pushViewport(next.vp)
      grid.yaxis(at= seq(0, 1, length= length(lab)), label= lab)
    popViewport()
	
    ## Plot X axis labels.
    next.vp <- viewport(x= 0.1, y= 0, width = 0.9, height = 0.15, just=c("left","bottom"))
    pushViewport(next.vp)
    for(i in 1:NROW(signal)) {
      grid.text(label[i], x= centers[i], y= 0.95, rot = 65, just=c("right", "top"))
    }
    popViewport()

  }

drawBarsVertical <- function(signal, error, label) {
    ## Remove NAs.
    indx <- !is.na(signal)
    signal <- signal[indx]
    error <- error[indx]
    label <- label[indx]

    ## compute other parts...
    lab<- pretty(signal, n=5)
    yrange <- max(lab, na.rm=TRUE) - min(lab, na.rm=TRUE)
    width <- 1/NROW(signal)/3 
    centers <- seq(width, 1, 1/NROW(signal))
    zero <- 1-max(lab, na.rm=TRUE)/yrange

    ord <- order(signal)
    signal <- signal[ord]
    error  <- error[ord]
    label  <- label[ord]
	
    ## Print the barplot.
    require(grid)
    grid.newpage()
    top.vp <- viewport(width = 0.98, height = 0.98, xscale= c(0,1), yscale= c(0,1))
    pushViewport(top.vp)
    next.vp <- viewport(y= 0.1, x= 0.15, height = 0.9, width = 0.85, just=c("left","bottom"))
    pushViewport(next.vp)
    for(i in lab) {
      grid.lines(x= (i-min(lab))/yrange, gp=gpar(col="light gray"))
    }
  
    for(i in 1:NROW(signal)) {
      x <- c(centers[i]+width, centers[i]-width, centers[i]-width, centers[i]+width)
      ym <- signal[i]-min(lab)
      y <- c(zero, zero, c(ym, ym)/yrange)
	  
      if(zero>(ym/yrange)) {
         fill="light blue"
         ye <- ym-error[i]
      }
      else {
         fill="dark red"
         ye <- ym+error[i]
      }

      grid.polygon(y, x, gp=gpar(fill=fill))
	
   ## Draw std. error.
      grid.lines(y=centers[i], x=c(ym, ye)/yrange)
      grid.lines(y=x[1:2], x=(ye)/yrange)
    }
    popViewport()

    ## Plot Y axis.
    next.vp <- viewport(y= 0.07, x= 0.15, height = 0.03, width = 0.85, just=c("left","bottom"))
    pushViewport(next.vp)
      grid.xaxis(at= seq(0, 1, length= length(lab)), label= lab)
    popViewport()
	
    ## Plot X axis labels.
    next.vp <- viewport(y= 0.1, x= 0, height = 0.9, width = 0.15, just=c("left","bottom"))
    pushViewport(next.vp)
    for(i in 1:NROW(signal)) {
      grid.text(label[i], y= centers[i], x= 0.95, rot = 0, just=c("right", "top"))
    }
    popViewport()

  }

  

