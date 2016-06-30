library(dREG)

gdm=genomic_data_model(c(10, 25, 50, 500, 5000), c(10, 10, 30, 20, 20))

bw_plus  ='/local/storage/data/hg19/k562/proseq/K562_unt.sort.bed.gz_minus.bw'
bw_minus ='/local/storage/data/hg19/k562/proseq/K562_unt.sort.bed.gz_plus.bw'

bed_file <- "/local/storage/data/hg19/k562/dreg/k562.bg.gz"

bed <- read.table(pipe(paste("zcat ", bed_file)), sep="\t");

getTesting<-function(bed, bw_p, bw_m, batch_size=50000, ncores=3){
    test = read_genomic_data(gdm, bed[,c(1,2,3)], bw_p, bw_m, batch_size=batch_size, ncores=ncores)
    print(dim(test))
    return(test);
}

time.n12.b50<- c();
time.n32.b50<- c();
time.n32.b5 <- c();

time.n32.b50<- rbind(time.n32.b50, c(10000, 32, 50, c( system.time(test <- getTesting(bed[1:10000,], bw_plus,bw_minus, batch_size=50000, ncores=32)) ) ))

time.n32.b50<- rbind(time.n32.b50, c(100000, 32, 50, c( system.time(test <- getTesting(bed[1:100000,], bw_plus,bw_minus, batch_size=50000, ncores=32)) )))

time.n32.b50<- rbind(time.n32.b50, c(1000000, 32, 50, c( system.time(test <- getTesting(bed[1:1000000,], bw_plus,bw_minus, batch_size=50000, ncores=32)) )))

time.n32.b50<- rbind(time.n32.b50, c(10000000, 32, 50, c( system.time(test <- getTesting(bed[1:10000000,], bw_plus,bw_minus, batch_size=50000, ncores=32)) ))); rm(test); gc();

rm(test); gc();

save(time.n32.b5, time.n12.b50, time.n32.b50,  file='benchmark-test-2016.rdata')

time.n32.b5 <- c();

time.n32.b5<- rbind(time.n32.b5, c(10000, 32, 5, c( system.time(test <- getTesting(bed[1:10000,], bw_plus,bw_minus, batch_size=5000, ncores=32)) )))

time.n32.b5<- rbind(time.n32.b5, c(100000, 32, 5, c( system.time(test <- getTesting(bed[1:100000,], bw_plus,bw_minus, batch_size=5000, ncores=32)) )))

time.n32.b5<- rbind(time.n32.b5, c(1000000, 32, 5, c( system.time(test <- getTesting(bed[1:1000000,], bw_plus,bw_minus, batch_size=5000, ncores=32)) )))

time.n32.b5<- rbind(time.n32.b5, c(10000000, 32, 5, c( system.time(test <- getTesting(bed[1:10000000,], bw_plus,bw_minus, batch_size=5000, ncores=32)) )))
rm(test); gc();

save(time.n32.b5, time.n12.b50, time.n32.b50, file='benchmark-test-2016.rdata')

time.n12.b50<- c();

time.n12.b50<- rbind(time.n12.b50, c(10000, 12, 50, c( system.time(test <- getTesting(bed[1:10000,], bw_plus,bw_minus, batch_size=50000, ncores=12)) )))

time.n12.b50<- rbind(time.n12.b50, c(100000, 12, 50, c( system.time(test <- getTesting(bed[1:100000,], bw_plus,bw_minus, batch_size=50000, ncores=12)) )))

time.n12.b50<- rbind(time.n12.b50, c(1000000, 12, 50, c( system.time(test <- getTesting(bed[1:1000000,], bw_plus,bw_minus, batch_size=50000, ncores=12)) )))

time.n12.b50<- rbind(time.n12.b50, c(10000000, 12, 50, c( system.time(test <- getTesting(bed[1:10000000,], bw_plus,bw_minus, batch_size=50000, ncores=12)) )))
rm(test); gc();

save(time.n32.b5, time.n12.b50, time.n32.b50, file='benchmark-test-2016.rdata')


time.n32.b500 <- c();
time.n32.b500<- rbind(time.n32.b500, c(10000, 32, 500, c( system.time(test <- getTesting(bed[1:10000,], bw_plus,bw_minus, batch_size=500000, ncores=32)) ) ))
time.n32.b500<- rbind(time.n32.b500, c(100000, 32, 500, c( system.time(test <- getTesting(bed[1:100000,], bw_plus,bw_minus, batch_size=500000, ncores=32)) )))
time.n32.b500<- rbind(time.n32.b500, c(1000000, 32, 500, c( system.time(test <- getTesting(bed[1:1000000,], bw_plus,bw_minus, batch_size=500000, ncores=32)) )))
time.n32.b500<- rbind(time.n32.b500, c(10000000, 32, 500, c( system.time(test <- getTesting(bed[1:10000000,], bw_plus,bw_minus, batch_size=500000, ncores=32)) )));
rm(test); gc();

save(time.n12.b50, time.n32.b5, time.n32.b50, time.n32.b500, file='benchmark-test-2016.rdata')

perform_fig<-function()
{
	load("benchmark-test-2016.rdata")
	t50.12 <- time.n12.b50[,1]/(time.n12.b50[,6]*time.n12.b50[,2]);
	t5 <- time.n32.b5[,1]/(time.n32.b5[,6]*time.n32.b5[,2]);
	t50 <- time.n32.b50[,1]/(time.n32.b50[,6]*time.n32.b50[,2]);
	t500 <- time.n32.b500[,1]/(time.n32.b500[,6]*time.n32.b500[,2]);
	x <- log10(time.n32.b50[,1]);
	y.max <- max(c(t5, t50, t500, t50.12))*1.25;

	load("benchmark-test-2015.rdata")
	tt50.12 <- time.n12.b50[,1]/(time.n12.b50[,6]*time.n12.b50[,2]);
	tt5 <- time.n32.b5[,1]/(time.n32.b5[,6]*time.n32.b5[,2]);
	tt50 <- time.n32.b50[,1]/(time.n32.b50[,6]*time.n32.b50[,2]);

	plot(x, t500, type="l", ylim=c(0, y.max),col="red", xlab="Genomic loci", ylab="Loci/second/cpu", xaxt="n", lwd=2);
	axis(1,at=c(4,5,6,7), labels=c("10^4", "10^5", "10^6", "10^7") )
	lines(x, t50, col="green", lwd=2);
	lines(x, t5, col="blue", lwd=2);
	lines(x, t50.12, col="purple", lwd=0.5);

	lines(x, tt5, col="blue", lty="22", lwd=2);
	lines(x, tt50, col="green", lty="22", lwd=2);
	lines(x, tt50.12, col="purple", lty="22", lwd=1/2);

	legend("topleft",
		legend=c("Batchsize=500/32 Cores/2016", "Batchsize=50/32 Cores/2016", "Batchsize=5/32 Cores/2016", "Batchsize=50/12 Cores/2016", "Batchsize=50/32 Cores/2015", "Batchsize=5/32 Cores/2015","Batchsize=50/12 Cores/2015"),
		col=c("red", "green", "blue", "purple", "green", "blue", "purple"),
		lwd=c(2,2,2,1/2,2,2,1/2),
		lty=c("solid","solid","solid","solid", "22","22","22"));

	dev.off();
}
