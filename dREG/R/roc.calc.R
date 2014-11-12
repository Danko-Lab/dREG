#
#  Creates a ROC plot, given: 
#    (1) A set of 'true' genomic intervals (e.g. ChIP-seq peaks).
#    (2) A set of 'possible' genomic intervals (e.g. DNAse-1 peaks).
#    (3) A vector of scores for each possibe genomic interval in (#2).
#

roc.calc <- function(true, possible, scores, filterPossible=TRUE, n_points= 100) {

## Get subset of 'true' peaks that overlap at least 1 'possible' peak.
 true.feat <- feat(seqname= true[,1], start= true[,2], end= true[,3])
 possible.feat <- feat(seqname= possible[,1], start= possible[,2], end= possible[,3])
 ol <- overlap.feat(x= true.feat, filter= possible.feat)
 indx <- match(paste(ol$seqname, ol$start, ol$end), paste(true.feat$seqname, true.feat$start, true.feat$end))
 if(filterPossible) {
   true.feat <- true.feat[indx,]
 }
 
## Get a good set of threshold values.
 prune.threshold <- n_points
 vunq <- unique(sort(scores))
 if(NROW(vunq)>prune.threshold) # Prune it down a bit...
	vunq <- vunq[c(round(seq(2, NROW(vunq), NROW(vunq)/(prune.threshold-2))), (NROW(vunq)-1))] ## Can't do 1 or NROW(vunq), becuase these can be -Inf or Inf.
 
## Foreach score in vunq, get TP, TN, FP, and FN.
TPR <- double()
FPR <- double()
 for(i in c(1:NROW(vunq))) {
  pos.indx <- scores > vunq[i]
  
   ## TP (true positives).  Number pos. that overlap true.feat.
   pos.feat <- feat(seqname= possible[pos.indx,1], start= possible[pos.indx,2], end= possible[pos.indx,3])
   TP <- NROW(overlap.feat(x= true.feat, filter= pos.feat))
   ## FN (false negatives).  #Num true features that I call negative ...
   FN <- NROW(true.feat)-TP #NROW(overlap.feat(x= neg.feat, filter= true.feat))

   ## These calculations are a bit unfair ... when filterPossible==FALSE ...
   #    I'm using the entire genome to  to compute the sensitivity.  
   #    But only those portions that I call to compute specificity.
   #    In theory ... I should be using most of the genome for all computations ...
   
   ## FP (false positives).  #pos w/ no overlap.
   FP <- NROW(overlap.feat(x= pos.feat, filter= true.feat, overlapping=FALSE))
   ## TN (true negatives).  Number !pos. that !overlap true.feat.
   neg.feat <- feat(seqname= possible[!pos.indx,1], start= possible[!pos.indx,2], end= possible[!pos.indx,3])
   TN <- NROW(overlap.feat(x= neg.feat, filter= true.feat, overlapping=FALSE)) 
   
   
   TPR[i] <- TP/(TP+FN)
   FPR[i] <- FP/(FP+TN)
 }

 return(data.frame(FPR= FPR, TPR= TPR, threshold= vunq))
}

## true -- a boolean
## scores -- thresholded to get a boolean
logreg.roc.calc <- function(true, scores) {
 ## Get a good set of threshold values.
 prune.threshold <- 100
 vunq <- unique(sort(scores))
 if(NROW(vunq)>prune.threshold) # Prune it down a bit...
	vunq <- vunq[c(round(seq(2, NROW(vunq), NROW(vunq)/(prune.threshold-2))), (NROW(vunq)-1))] ## Can't do 1 or NROW(vunq), becuase these can be -Inf or Inf.
 
 ## Foreach score in vunq, get TP, TN, FP, and FN.
 TPR <- double()
 FPR <- double()
 for(i in c(1:NROW(vunq))) {
  pos.indx <- scores > vunq[i]
  
   ## TP (true positives).  Number pos. that overlap true.feat.
   TP <- sum(true & pos.indx)
   
   ## FP (false positives).  #pos w/ no overlap.
   FP <- sum(!true & pos.indx)

   ## FN (false negatives).  Number !pos. that overlap true.feat.
   FN <- sum(true & !pos.indx)
   
   ## TN (true negatives).  #NOT pos w/ no overlap.
   TN <- sum(!true & !pos.indx)
   
   TPR[i] <- TP/(TP+FN)
   FPR[i] <- FP/(FP+TN)
 }

 return(data.frame(FPR= FPR, TPR= TPR, threshold= vunq))
}

## Combines ROC plots, interpolating and weighting by nTP.
combine.roc <- function(list.roc, weight=rep(1, NROW(list.roc)), interp.corners=FALSE, use.max=FALSE, nvals=100) {
  FPR <- seq(0, 1, 1/nvals)
  
  ## Interpolate each 
  for(i in 1:NROW(list.roc)) {
    if(interp.corners) list.roc[[i]] <- rbind(c(1, 1, 1), list.roc[[i]], c(0, 0, 0))
	interp.tf <- sapply(FPR, function(x) {
	  ## Find closest two points.
	  indxMin <- min(which(list.roc[[i]]$FPR < x))
	  indxMax <- max(which(list.roc[[i]]$FPR > x))
	  
	  ## Handle out of defined region.
	  if(is.infinite(indxMin) | is.infinite(indxMax)) return(NA)
	  
	  ## Interpolate.
	  m <- (list.roc[[i]]$TPR[indxMax]-list.roc[[i]]$TPR[indxMin])/(list.roc[[i]]$FPR[indxMax]-list.roc[[i]]$FPR[indxMin])
	  b <- (list.roc[[i]]$TPR[indxMin])-m*list.roc[[i]]$FPR[indxMin]
	  
	  ## Return FPR.
	  return(m*x+b)
	})
	if(i == 1) interpProd <- interp.tf
	else interpProd <- cbind(interpProd, interp.tf)
  }

  TPR <- sapply(c(1:NROW(FPR)), function(x) { 
	if(use.max) {w.TPR <- max(interpProd[x,], na.rm=TRUE)
	} else {w.TPR <- sum(interpProd[x,]*weight, na.rm=TRUE)/sum(weight[!is.na(interpProd[x,])]) }
	if(is.nan(w.TPR)) {return(NA)
	} else { if(is.infinite(w.TPR)) {return(NA)
	} else { return(w.TPR) }}
  } )
  ## Multiply vectors, weighting 
  
  return(data.frame(FPR, TPR)[!is.na(TPR),])
}

## Computes the AUC of a ROC plot.
# TEST: 
#roc.auc( data.frame(FPR=c(0, 0.25, 0.5, 0.75, 1), TPR=c(0, 0.5, 0.8, 0.95, 1), threshold=c(1, 1, 1, 1, 1)))
#roc.auc( data.frame(FPR=c(0, 0.25, 0.5, 0.75, 1), TPR=c(0, 0.25, 0.5, 0.75, 1), threshold=c(1, 1, 1, 1, 1)))
#roc.auc( data.frame(FPR=c(0, 0.25, 0.5, 0.75, 1), TPR=c(0, 0.55, 0.85, 0.95, 1), threshold=c(1, 1, 1, 1, 1)))
roc.auc <- function(ROC) {
  ROC <- ROC[order(ROC$FPR),]
  area <- 0.5*(0+ROC$TPR[1])*(ROC$FPR[1]-0)
  for(i in 1:(NROW(ROC)-1)) {
    area <- area+ 0.5*(ROC$TPR[i]+ROC$TPR[i+1])*(ROC$FPR[i+1]-ROC$FPR[i])## Trapezoid area.
  }
  area <- area+ 0.5*(ROC$TPR[NROW(ROC)]+1)*(1-ROC$FPR[NROW(ROC)])
  return(area)
}

roc.plot <- function(ROC, ...) {
  plot(ROC[,c(1:2)], xlab="False positive rate", ylab="True positive rate", type="l", ...) #
  abline(0, 1, lty="dashed", col="gray", xlim=c(0,1), ylim=c(0,1))
}


