#make metaplots for S and U from GRO-cap labels




## Read in bigWigs.
require(bigWig)


# modify plot.metaprofile function to overlay 4 lines on each plot

plot.metaprofile.multiple<-function (x.list, minus.profile = NULL, X0 = plus.profile$X0, draw.error = TRUE, 
    cols = c("red", "blue", "lightgrey", "lightgrey"), ylim = NULL, 
    xlim = NULL, xlab = "Distance (bp)", ylab = plus.profile$name, alpha=0.3,
    ...) 
{
		
		
		stopifnot(length(x.list)==length(cols))
	
	  col = c(cols[1], "blue", "lightgrey")
	  
	  ploygon.col<- rgb (red= col2rgb(col[3])[1]/255,green= col2rgb(col[3])[2]/255,blue=col2rgb(col[3])[3]/255,alpha= alpha)
	  
	  
	  	  
	  ymin_vec<-c()
	  for (i in 1:length(x.list)){
	  	ymin_vec<-c(ymin_vec, x.list[[i]]$bottom)
	  }
	  
	   ymax_vec<-c()
	  for (i in 1:length(x.list)){
	  	ymax_vec<-c(ymax_vec, x.list[[i]]$top)
	  }
	  

	  
	  ymax<-max(ymax_vec)
	  ymin<-min(ymin_vec)
	  
	  
	  ylim=c(ymin,ymax)
	  
	  	  
		x<-x.list[[1]]
    plus.profile = x
    opt.args = list(...)
    if (!(all(plus.profile$top >= 0) && all(plus.profile$middle >= 
        0) && all(plus.profile$bottom >= 0))) 
        warning("expected plus profile curves to all be >= 0")
    if (!is.null(minus.profile) && !(all(minus.profile$top >= 
        0) && all(minus.profile$middle >= 0) && all(minus.profile$bottom >= 
        0))) 
        warning("expected minus profile curves to all be >= 0")
    N = length(plus.profile$middle)
    step = plus.profile$step
    if (is.null(step)) 
        stop("plus.profile is missing the 'step' element")
    stopifnot(X0 <= length(plus.profile$middle) * step && X0 >= 
        0)
    if (!is.null(minus.profile)) {
        stopifnot(length(plus.profile$middle) == length(minus.profile$middle))
        stopifnot(plus.profile$step == minus.profile$step)
    }
    x = 1:N * step - X0
    if (is.null(xlim)) 
        xlim = c(min(x), max(x))
    if (is.null(ylim)) {
        if (draw.error) {
            if (!is.null(minus.profile)) {
                ylim = c(-max(minus.profile$top, minus.profile$middle, 
                  minus.profile$bottom), max(plus.profile$top, 
                  plus.profile$middle, plus.profile$bottom))
            }
            else {
                ylim = c(min(plus.profile$top, plus.profile$middle, 
                  plus.profile$bottom), max(plus.profile$top, 
                  plus.profile$middle, plus.profile$bottom))
            }
        }
        else {
            if (!is.null(minus.profile)) {
                ylim = c(-max(minus.profile$middle), max(plus.profile$middle))
            }
            else {
                ylim = c(min(plus.profile$middle), max(plus.profile$middle))
            }
        }
    }
    plot(x, plus.profile$middle, col = col[1], xlim = xlim, ylim = ylim, 
        type = "l", xlab = xlab, ylab = ylab, ...)
    if (draw.error) {
        polygon(c(x, rev(x)), c(plus.profile$top, rev(plus.profile$bottom)), 
            col = ploygon.col, border = NA, ...)
        if (!is.null(minus.profile)) 
            polygon(c(x, rev(x)), -c(minus.profile$bottom, rev(minus.profile$top)), 
                col = ploygon.col, border = NA, ...)
    }
    lines(x, plus.profile$middle, col = col[1], ...)
    if (!is.null(minus.profile)) {
        lines(x, -minus.profile$middle, col = col[2], ...)
        lines(x, rep(0, N), ...)
    }
    
    
        
    for(i in c(2:length(x.list))){
    	
     col = c(cols[i], "blue", "lightgrey", "lightgrey")

    	x<-x.list[[i]]
    	plus.profile = x
    opt.args = list(...)
    if (!(all(plus.profile$top >= 0) && all(plus.profile$middle >= 
        0) && all(plus.profile$bottom >= 0))) 
        warning("expected plus profile curves to all be >= 0")
    if (!is.null(minus.profile) && !(all(minus.profile$top >= 
        0) && all(minus.profile$middle >= 0) && all(minus.profile$bottom >= 
        0))) 
        warning("expected minus profile curves to all be >= 0")
    N = length(plus.profile$middle)
    step = plus.profile$step
    if (is.null(step)) 
        stop("plus.profile is missing the 'step' element")
    stopifnot(X0 <= length(plus.profile$middle) * step && X0 >= 
        0)
    if (!is.null(minus.profile)) {
        stopifnot(length(plus.profile$middle) == length(minus.profile$middle))
        stopifnot(plus.profile$step == minus.profile$step)
    }
    x = 1:N * step - X0
    if (is.null(xlim)) 
        xlim = c(min(x), max(x))
    if (is.null(ylim)) {
        if (draw.error) {
            if (!is.null(minus.profile)) {
                ylim = c(-max(minus.profile$top, minus.profile$middle, 
                  minus.profile$bottom), max(plus.profile$top, 
                  plus.profile$middle, plus.profile$bottom))
            }
            else {
                ylim = c(min(plus.profile$top, plus.profile$middle, 
                  plus.profile$bottom), max(plus.profile$top, 
                  plus.profile$middle, plus.profile$bottom))
            }
        }
        else {
            if (!is.null(minus.profile)) {
                ylim = c(-max(minus.profile$middle), max(plus.profile$middle))
            }
            else {
                ylim = c(min(plus.profile$middle), max(plus.profile$middle))
            }
        }
    }
    lines(x, plus.profile$middle, col = col[1], xlim = xlim, ylim = ylim, 
        type = "l", xlab = xlab, ylab = ylab, ...)
    if (draw.error) {
        polygon(c(x, rev(x)), c(plus.profile$top, rev(plus.profile$bottom)), 
            col = ploygon.col, border = NA, alpha= alpha, ...)
        if (!is.null(minus.profile)) 
            polygon(c(x, rev(x)), -c(minus.profile$bottom, rev(minus.profile$top)), 
                col = ploygon.col, border = NA, alpha= alpha, ...)
    }
    lines(x, plus.profile$middle, col = col[1], ...)
    if (!is.null(minus.profile)) {
        lines(x, -minus.profile$middle, col = col[2], ...)
        lines(x, rep(0, N), ...)
    }
    }

}




doit_SU<-function(S.plus.bed,S.minus.bed,U.plus.bed,U.minus.bed,stp, coverage,target.bw,cols=c("red","blue"), ...) {
	
	scaling_func<-function(x){
		x/ stp
	}
	
	
	
	S.plus.bed <-read.table(S.plus.bed)
	S.minus.bed <-read.table(S.minus.bed)
	U.plus.bed <-read.table(U.plus.bed)
	U.minus.bed <-read.table(U.minus.bed)
	
	
	S.plus.bed<-center.bed(S.plus.bed, (coverage-1), coverage)
	U.plus.bed<-center.bed(U.plus.bed, (coverage-1), coverage)
	S.minus.bed <-center.bed(S.minus.bed, coverage, (coverage-1))	
	U.minus.bed <-center.bed(U.minus.bed, coverage, (coverage-1))	
	
	S.bed<-rbind.data.frame(S.plus.bed, S.minus.bed)
	U.bed<-rbind.data.frame(U.plus.bed, U.minus.bed)
	
	
	
	target.bw<-load.bigWig(target.bw)
	
	
	S.profile <- metaprofile.bigWig(S.bed, bw.plus=target.bw, bw.minus=target.bw, step=stp, matrix.op= scaling_func)
	U.profile <- metaprofile.bigWig(U.bed, bw.plus=target.bw, bw.minus=target.bw, step=stp, matrix.op= scaling_func)
	
	
	x.list<-list(S.profile, U.profile)
	
	unload.bigWig(target.bw)

	plot.metaprofile.multiple(x.list= x.list,X0= coverage/stp, cols = cols,... )

        
}



if(0)
{


coverage <- 5000

stp=50


markers<-list.files(path="/workdir/tc532/sep_classifier/k562/histones",pattern="*.bigWig")

markers_path<-list.files(path="/workdir/tc532/sep_classifier/k562/histones",pattern="*.bigWig",full.names=T)

markers_name <-as.character(sapply(markers_path,FUN=function(x) unlist(strsplit (unlist(strsplit(x,split="K562"))[2],split="StdSig"))[1]))


bed.path<-"/workdir/tc532/sep_classifier/k562/training.bed"
#bed.path<-"training.bed"

S.plus.bed=paste(bed.path,"hg19.k562.new_hmm2b.post2.S_plus.bdRM.bed",sep="/")
S.minus.bed=paste(bed.path,"hg19.k562.new_hmm2b.post2.S_minus.bdRM.bed",sep="/")
U.plus.bed=paste(bed.path,"hg19.k562.new_hmm2b.post2.U_plus.bdRM.bed",sep="/")
U.minus.bed=paste(bed.path,"hg19.k562.new_hmm2b.post2.U_minus.bdRM.bed",sep="/")



pdf(paste("metaplot.", coverage, ".pdf",sep=""),width=12.444, height=7,pointsize=5,useDingbats=FALSE )
for (k in 1:length(markers)){
	print(markers_name[k])
	
	doit_SU(S.plus.bed,S.minus.bed,U.plus.bed,U.minus.bed,stp,coverage,target.bw=markers_path[k],cols=c("red","blue"),main= markers_name[k])
}

dev.off()


}




if(0)
{




coverage <- 50000

stp=50


markers<-list.files(path="/workdir/tc532/sep_classifier/gm12878/histones",pattern="*.bigWig")

markers_path<-list.files(path="/workdir/tc532/sep_classifier/gm12878/histones",pattern="*.bigWig",full.names=T)

markers_name <-as.character(sapply(markers_path,FUN=function(x) unlist(strsplit (unlist(strsplit(x,split="gm12878"))[2],split="StdSig"))[1]))


bed.path<-"/workdir/tc532/sep_classifier/gm12878/training.bed"
#bed.path<-"training.bed"

S.plus.bed=paste(bed.path,"hg19.gm12878.new_hmm2b.post2.S_plus.bdRM.bed",sep="/")
S.minus.bed=paste(bed.path,"hg19.gm12878.new_hmm2b.post2.S_minus.bdRM.bed",sep="/")
U.plus.bed=paste(bed.path,"hg19.gm12878.new_hmm2b.post2.U_plus.bdRM.bed",sep="/")
U.minus.bed=paste(bed.path,"hg19.gm12878.new_hmm2b.post2.U_minus.bdRM.bed",sep="/")



pdf(paste("metaplot.gm12878.", coverage, ".pdf",sep=""),width=12.444, height=7,pointsize=5,useDingbats=FALSE )
for (k in 1:length(markers)){
	print(markers_name[k])
	
	doit_SU(S.plus.bed,S.minus.bed,U.plus.bed,U.minus.bed,stp,coverage,target.bw=markers_path[k],cols=c("red","blue"),main= markers_name[k])
}

dev.off()



}





