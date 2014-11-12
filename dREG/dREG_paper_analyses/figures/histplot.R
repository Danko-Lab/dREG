## Lollypop plot type thingamig.

cd.circle <- function(x, y, r) {
  rads <- pi*c(32:0)/64
  qdx <- r*cos(rads)
  qdy <- r*sin(rads)
  dx <- c(qdx, rev(qdx), -1*qdx, rev(-1*qdx))
  dy <- c(qdy, rev(-1*qdy), -1*qdy, rev(qdy))
  return(list(x= x+dx, y= y+dy))
} ## TEST!
# cd.circle(0.5, 0.5, 0.5)
# require(grid)
# grid.newpage()
# top.vp <- viewport(width = 0.98, height = 0.98, xscale= c(0,1), yscale= c(0,1))
# pushViewport(top.vp)
# cc <- cd.circle(0.5, 0.5, 0.5)
# grid.polygon(cc$x, cc$y)

histplot<- function(histones, fillCols) {
  lab <- pretty(c(0,1), n=5)

  ## Reorder.
  histones <- histones[order(sapply(histones, function(x) {mean(x[[2]])}))]
  
  ## Newpage.
  require(grid)
  grid.newpage()
  top.vp <- viewport(width = 0.98, height = 0.98, xscale= c(0,1), yscale= c(0,1))
  pushViewport(top.vp)

  ## Plot histone bars/ ... .
  width <- 1/NROW(histones)/3
  centers <- seq(width, 1, 1/NROW(histones))

  next.vp <- viewport(x= 0.1, y= 0.15, width = 0.9, height = 0.85, just=c("left","bottom"))
  pushViewport(next.vp)
  
  ## Gridlines.
  for(i in lab) {
    grid.lines(y= i, gp=gpar(col="light gray"))
  }
  
  for(i in 1:NROW(histones)) {
    x <- c(centers[i]+width, centers[i]-width, centers[i]-width, centers[i]+width)
	ym <- mean(histones[[i]][[2]])/100
    y <- c(0, 0, ym, ym)
 #   grid.polygon(x, y, gp=gpar(fill=histones[[i]][[4]]))
    grid.lines(centers[i], c(0, ym))
	npoints <- NROW(histones[[i]][[2]])
	for(j in 1:npoints) {
	  cc <- cd.circle(centers[i], histones[[i]][[2]][j]/100, width/2)
	  grid.polygon(cc$x+(width/4)-(width*j/npoints/4), cc$y, gp=gpar(fill=fillCols[[names(histones[[i]][[2]])[j]]]))
	}
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
  for(i in 1:NROW(histones)) {
    grid.text(histones[[i]][[1]], x= centers[i], y= 0.95, rot = 65, just=c("right", "top"))
  }
  popViewport()

}
#histplot(histones)


cd.lolyplot<- function(data, names, fill="black") {
  lab <- pretty(c(0,1), n=5)
  labels <- unique(names)
  
  ## Newpage.
  require(grid)
  grid.newpage()
  top.vp <- viewport(width = 0.98, height = 0.98, xscale= c(0,1), yscale= c(0,1))
  pushViewport(top.vp)

  ## Plot histone bars/ ... .
  width <- 1/NROW(labels)/20
  centerseq <- seq(width, 1, 1/NROW(labels))
  centers <- centerseq[as.integer(as.factor(names))]
  
  next.vp <- viewport(x= 0.1, y= 0.15, width = 0.9, height = 0.85, just=c("left","bottom"))
  pushViewport(next.vp)
  
  ## Gridlines.
  for(i in lab) {
    grid.lines(y= i, gp=gpar(col="light gray"))
  }
  
  for(i in 1:NROW(data)) {
	ym <- data[i]
	cc1 <- cd.circle(centers[i], ym, width)
	grid.polygon(cc1$x, cc1$y, gp=gpar(fill=fill))
  }
  
  ## Medians/ Means
  for(i in 1:NROW(labels)) {
    indx <- levels(as.factor(names))==labels[i]
    grid.lines(centerseq[i]+c(-1/12/NROW(labels),1/12/NROW(labels)), mean(data[names==labels[indx]])) ## Draw horizontal line at mean.
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
  for(i in 1:NROW(labels)) {
    grid.text(labels[levels(as.factor(names))==labels[i]], x= centerseq[i], y= 0.95, rot = 65, just=c("right", "top"))
  }
  popViewport()

}
## test: 
## cd.lolyplot(c(1:10)/10, c(rep("A",4), rep("B", 2), rep("C", 4)))

cd.barplot<- function(data, error, names, fill) {
  lab <- pretty(c(0,1), n=5)

  ## Reorder.
  ord <- order(data)
  data <- data[ord]
  names <- names[ord]
  error <- error[ord]
  
  ## Newpage.
  require(grid)
  grid.newpage()
  top.vp <- viewport(width = 0.98, height = 0.98, xscale= c(0,1), yscale= c(0,1))
  pushViewport(top.vp)

  ## Plot histone bars/ ... .
  width <- 1/NROW(data)/3
  centers <- seq(width, 1, 1/NROW(data))

  next.vp <- viewport(x= 0.1, y= 0.15, width = 0.9, height = 0.85, just=c("left","bottom"))
  pushViewport(next.vp)
  
  ## Gridlines.
  for(i in lab) {
    grid.lines(y= i, gp=gpar(col="light gray"))
  }
  
  for(i in 1:NROW(data)) {
	ym <- data[i]
    x <- c(centers[i]+width, centers[i]-width, centers[i]-width, centers[i]+width)
    y <- c(0, 0, ym, ym)
    grid.polygon(x, y, gp=gpar(fill=fill))

   ## Draw std. error.
    grid.lines(x=centers[i], y=c(ym, ym+error[i]))
    grid.lines(x=x[1:2], y=(ym+error[i]))
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
  for(i in 1:NROW(data)) {
    grid.text(names[i], x= centers[i], y= 0.95, rot = 65, just=c("right", "top"))
  }
  popViewport()

}
