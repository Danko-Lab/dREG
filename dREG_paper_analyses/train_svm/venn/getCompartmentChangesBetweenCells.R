##
## Evaluates changes between compartments, and compares across cell types.
##

createTmp <- function(f, tmpName) {
  fnames= c("chromHMM.only.bed", "plus.DNAse.bed", "plus.dREG.bed", "plus.dREG.PROONLY.bed", "plus.dREG.ENHONLY.bed", "insulator.bed")
  sizes <- integer(0)
  tab <- NULL
  
  ## Read all elements and append to a long list.
  for(i in fnames) {
    tmp <- read.table(paste(f,i,sep="."))[,1:3]
    sizes <- c(sizes,  NROW(tmp))
	tab <- rbind(tab, tmp)
  }
  
  ## Write each file.
  for(i in 1:NROW(fnames)) {
    samp <- sample(c(1:NROW(tab)), sizes[i], replace=FALSE)
	revSamp <- which(!(1:NROW(tab) %in% samp))
    out_tab <- tab[samp,]  ## Create a sample to write out.
	tab <- tab[revSamp,]   ## Remove those from the dataset.
	
	write.table(out_tab, pipe(paste("sort-bed - > ",tmpName,".",fnames[i], sep="")), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
  }
}

## a,b -> k562|gm12878|hela|cd4
labelswap <- function(a, b, n=2) {
  MAT <- list()
  for(i in 1:n) {
    #createTmp(a, "tmp1")
	createTmp(b, "tmp2")
    MAT[[i]] <- read.table(text= system(paste("bash getCompartmentChangesBetweenCells.bsh",a,"tmp2"), intern=TRUE))
	system("rm tmp2*") # tmp1* 
  }
  MAT
}

## Returns means of elements across a list of matrices/ data.frames (listObj).
listMeans <- function(listObj) {
  Reduce("+", listObj) / length(listObj)
}
## Returns variences of elements across a list of matrices/ data.frames (listObj).
listVars  <- function(listObj) {
 means <- listMeans(listObj)
 listMeans(lapply(listObj, function(x) {(x-means)^2})) ## Varience.
}

## a,b -> k562|gm12878|hela|cd4
read.dataset <- function(a, b) {
 rd <- read.table(text= system(paste("bash getCompartmentChangesBetweenCells.bsh", a, b), intern=TRUE))
 
 return(rd/rowSums(rd))
 
 ## This returns permuted data.
# permutList <- labelswap(a, b, n=100) ## NO SWAP.
# means <- listMeans(permutList)
# vars  <- listVars(permutList)
# 
# print(rd)
# (rd-means)/sqrt(vars)
}

 kg <- read.dataset("k562", "gm12878") #read.table(text= system("bash getCompartmentChangesBetweenCells.bsh k562 gm12878", intern=TRUE))#; kg <- kg/rowSums(kg)
 kc <- read.dataset("k562", "cd4") #kc <- read.table(text= system("bash getCompartmentChangesBetweenCells.bsh k562 cd4", intern=TRUE))#; kc <- kc/rowSums(kc)
 kh <- read.dataset("k562", "hela") #read.table(text= system("bash getCompartmentChangesBetweenCells.bsh k562 hela", intern=TRUE))#; kh <- kh/rowSums(kh)

 gk <- read.dataset("gm12878", "k562") #read.table(text= system("bash getCompartmentChangesBetweenCells.bsh gm12878 k562", intern=TRUE))#; gk <- gk/rowSums(gk)
 gc <- read.dataset("gm12878", "cd4") #read.table(text= system("bash getCompartmentChangesBetweenCells.bsh gm12878 cd4", intern=TRUE))#; gc <- gc/rowSums(gc)
 gh <- read.dataset("gm12878", "hela") #read.table(text= system("bash getCompartmentChangesBetweenCells.bsh gm12878 hela", intern=TRUE))#; gh <- gh/rowSums(gh)

 cg <- read.dataset("cd4", "gm12878") #read.table(text= system("bash getCompartmentChangesBetweenCells.bsh cd4 gm12878", intern=TRUE))#; cg <- cg/rowSums(cg)
 ck <- read.dataset("cd4", "k562") #read.table(text= system("bash getCompartmentChangesBetweenCells.bsh cd4 k562", intern=TRUE))#; ck <- ck/rowSums(ck)
 ch <- read.dataset("cd4", "hela") #read.table(text= system("bash getCompartmentChangesBetweenCells.bsh cd4 hela", intern=TRUE))#; ch <- ch/rowSums(ch)

 hg <- read.dataset("hela", "gm12878") #read.table(text= system("bash getCompartmentChangesBetweenCells.bsh hela gm12878", intern=TRUE))#; hg <- hg/rowSums(hg)
 hk <- read.dataset("hela", "k562") #read.table(text= system("bash getCompartmentChangesBetweenCells.bsh hela k562", intern=TRUE))#; hk <- hk/rowSums(hk)
 hc <- read.dataset("hela", "cd4") #read.table(text= system("bash getCompartmentChangesBetweenCells.bsh hela cd4", intern=TRUE))#; hc <- hc/rowSums(hc)
 
 #save.image("cellSwitch/cellSwitch.RData")
 
 zm <- matrix(0, nrow=5, ncol=6); rownames(zm) <- rownames(kg); colnames(zm) <- colnames(kg)

 dafr <- cbind( k= rbind(zm,gm= gk,cd= ck),  gm= rbind(k= kg, zm, cd= cg), cd= rbind(k= kc, gm= gc, zm))
 
## This is a lot to wrap your head around. 
# levelplot(t(as.matrix(dafr))[c(NROW(t(dafr)):1),], col.regions=cols, scales = list(x= list(rot=90, y=list(draw=FALSE))))
 meanplot <- kg
 for(i in 1:NROW(kg))
   for(j in 1:NCOL(kg))
     meanplot[i,j] <- median(c(kg[i,j], kc[i,j], kh[i,j], gk[i,j], gc[i,j], gh[i,j], cg[i,j], ck[i,j], ch[i,j], hg[i,j], hk[i,j], hc[i,j])) #mean(kg[i,j], gk[i,j]) #
 meanplot
 #meanplot <- log(abs(meanplot))*(as.integer(meanplot>0)*2-1) ## Take a log, preserving the sign.

 pdf("cellSwitch/cellType.Heatmap.pdf")
  ## Heatmap ...  
  library(lattice)
  nbreaks=60
  cols= colorRampPalette(c("white", "#CA7779", "#DDCA3B"))(100) #c("#6600FF", "#FFFFFF", "#FF3333"))(nbreaks)# 
  peak= max(abs(meanplot))
  at= seq(0, peak, peak/nbreaks) #pretty(seq(-peak, peak), nbreaks)
  levelplot(t(as.matrix(meanplot)), col.regions=cols, scales = list(x= list(rot=90), y=list(draw=FALSE)), at=at)
 dev.off()
 
 TTp <- c(kg[4,5], kc[4,5], kh[4,5], gk[4,5], gc[4,5], gh[4,5], cg[4,5], ck[4,5], ch[4,5], hg[4,5], hk[4,5], hc[4,5])
 TTe <- c(kg[3,4], kc[3,4], kh[3,4], gk[3,4], gc[3,4], gh[3,4], cg[3,4], ck[3,4], ch[3,4], hg[3,4], hk[3,4], hc[3,4])
 DD <- c(kg[2,3], kc[2,3], kh[2,3], gk[2,3], gc[2,3], gh[2,3], cg[2,3], ck[2,3], ch[2,3], hg[2,3], hk[2,3], hc[2,3])
 CC <- c(kg[1,2], kc[1,2], kh[1,2], gk[1,2], gc[1,2], gh[1,2], cg[1,2], ck[1,2], ch[1,2], hg[1,2], hk[1,2], hc[1,2])
 
 wilcox.test(TTe, DD)
 wilcox.test(CC, DD)
 boxplot(TTp, TTe, DD, CC)
 
 Tpn <- c(kg[4,1], kc[4,1], kh[4,1], gk[4,1], gc[4,1], gh[4,1], cg[4,1], ck[4,1], ch[4,1], hg[4,1], hk[4,1], hc[4,1])
 Ten <- c(kg[3,1], kc[3,1], kh[3,1], gk[3,1], gc[3,1], gh[3,1], cg[3,1], ck[3,1], ch[3,1], hg[3,1], hk[3,1], hc[3,1])
 Dn <- c(kg[2,1], kc[2,1], kh[2,1], gk[2,1], gc[2,1], gh[2,1], cg[2,1], ck[2,1], ch[2,1], hg[2,1], hk[2,1], hc[2,1])
 Cn <- c(kg[1,1], kc[1,1], kh[1,1], gk[1,1], gc[1,1], gh[1,1], cg[1,1], ck[1,1], ch[1,1], hg[1,1], hk[1,1], hc[1,1])
 In <- c(kg[5,1], kc[5,1], kh[5,1], gk[5,1], gc[5,1], gh[5,1], cg[5,1], ck[5,1], ch[5,1], hg[5,1], hk[5,1], hc[5,1])

 wilcox.test(Ten, Dn)
 wilcox.test(Ten, Cn)
 wilcox.test(Cn, Dn)
 boxplot(Tpn, Ten, Dn, Cn)

 pdf("cellSwitch/cellTypeTransition.pdf")
  source("~/src/featureDetector/test_functions/figures/histplot.R")
#  cd.lolyplot( c(TT, DD, CC), c(rep("TT", NROW(TT)), rep("DD", NROW(DD)), rep("CC", NROW(CC))), "black")
  cd.lolyplot( c(Tpn, Ten, Dn, Cn, In), c(rep("Tpn", NROW(Tpn)), rep("Ten", NROW(Ten)), rep("Dn", NROW(Dn)), rep("Cn", NROW(Cn)), rep("In", NROW(In))), "black")
 dev.off()
 
