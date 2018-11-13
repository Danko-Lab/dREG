#                                                    P        S       E                           T   A  T  A                  
#RMRP		    	TTTTTTTAATCTCACGCCACCAACT -70 TTCTCACCCTAATCATAAAACA -49 CAATTTCTTTAG--GGC-33 TATAAAATACT -23 ACTCTGTGAAGCTGAGGACGTG
#RN7SK 		        TGCTGAAGCTCTAGTACGATAAGCA -68 ACTTGACC-TAAGTGTAAAGTT -48 GAGACTTCCTT---CAG-33 GTTTATATAGC -23 TTGTGCGCCGCTTGGGTACCTC
#RNU6 (U6-1)	    TTATGTTTTAAAATGGACTATCATA -68 TGCTTACCGTAACTTGAAAGTA -47 TTTCGATTTCT---TGG-32 CTTTATATATC -22 TTGTGGAAAGGACGAAACACC
#RNU6 (U6-2)	    TATGCTAAATATGAAACCGACCATA -67 AGTTATCC-TAACCAAAAGATG -47 ATTTGATTGA----AGG-33 GCTTAAAATAG -23 GTGTGACAGTAACCCTTGAGTC
#RNU6 (U6-7)	    AAGATGGACAGGAAAGCGCTCGATT -68 AGGTTACCGTAAGGAAAACAAA -47 TGAGAAACTCCC--GTG-31 CCTTATAAGAC -21 CTGGGGACGGACTTATTTGC
#RNU6 (U6-8)	    AAGATGGACAGGAAAGGGCGCGGTT -70 CGGTCACCGTAAGTAGAATAGG -49 TGAAAGACTCCC--GTG-33 CCTTATAAGGC -23 CTGTGGGTGACTTCTTCTCAAC
#RNU6 (U6-9)	    TTAAATCTCTAGGTCATTTAAGAGA -67 AGTCGGCC-TATGTGTACAGAC -47 ATTTGTTCCAG---GGG-32 CTTTAAATAGC -22 TGGTGGTGGAACTCAATATTC
#RNU6ATAC		    TGTACCTCCATGGATAGCGAACAAG -67 AAGTCACCCTCACCGAAAGGCG -46 AGTGGAGCTTT---CGT-31 CCTTAAATAAA -21 GTGCGCAGGGAAGCCGAGGC
#RNY1		    	TAGTCATCAGTAAACTGAAACCAGA -69 ATATCACTGTAAGGGGAAAATG -48 AACAAATTTGG---GGG-33 CTTTAAATAGT -23 TCAAACAGTAGGAGGACTTATT
#RNY3		    	CTTCCTTTTTTTTAGCTCCTGTGAA -70 TAGTCACCGTAACTATGGTAGA -49 GATGGAACTTTCGAGGC-31 TTATATAAGTA -21 GCAGCGTGCCTTTGTGTTTC
#RNY4			    GACTTTTTGGAGAATTCTTAAAATA -68 ACTCATCC-TAACTTATTTAGA -48 GTAGCCACTTCA--GAG-32 ATTTATAAAAT -22 GAAAGTGAAAGCAGTTTTTCT
#RNY5		    	TTACGATCATGGCATAGGCTCTGAA -67 AAGTCTCCTTACCTAGAAAAGA -46 CCCTAAGTAG----GCA-32 CTATAAATAAC -22 AAGAGACTCACAGGATAACAC 
#RPPH1		        TTTGCATGTCGCTATGTGTTCTGGG -71 AAATCACCATAAACGTGAAATG -50 TCTTTGGATTTGG-GAA-33 TCTTATAAGTT -23 CTGTATGAGACCACTTTTTCCC
#TRNAU1	   	        CAACCATCTCACACCTTTCCAAAGG -69 ACGCGACCATAACTCTAAAAGG -48 TAAGCTTTTGC---GAT-33 CCTTATATAGC -23 TGCGCGGGAATAAGGTTGTCCT


# PSE
PSE.seq <- c(
"TCACCCTAATCATAAAACA",
"TGACC-TAAGTGTAAAGTT",
"TTACCGTAACTTGAAAGTA",
"TATCC-TAACCAAAAGATG",
"TTACCGTAAGGAAAACAAA",
"TCACCGTAAGTAGAATAGG",
"CGGCC-TATGTGTACAGAC",
"TCACCCTCACCGAAAGGCG",
"TCACTGTAAGGGGAAAATG",
"TCACCGTAACTATGGTAGA",
"CATCC-TAACTTATTTAGA",
"TCTCCTTACCTAGAAAAGA",
"TCACCATAAACGTGAAATG")


#TATA
TATA.seq <- c("TATAAAATACT",
"GTTTATATAGC",
"CTTTATATATC",
"GCTTAAAATAG",
"CCTTATAAGAC",
"CCTTATAAGGC",
"CTTTAAATAGC",
"CCTTAAATAAA",
"CTTTAAATAGT",
"TTATATAAGTA",
"ATTTATAAAAT",
"CTATAAATAAC",
"TCTTATAAGTT",
"CCTTATATAGC")

library("tools")
options("scipen"=100, "digits"=4)
options(width=250)

bed_intersect <- function( score_bed, refer_bed ){
	file.refer <- tempfile(fileext=".bed");
	file.score <- tempfile(fileext=".bed");
	write.table( refer_bed[,c(1:3)], file=file.refer, quote=F, col.names=F, row.names=F, sep="\t");
	write.table( score_bed, file=file.score, quote=F, col.names=F, row.names=F, sep="\t");
	ret <- unique( read.table( pipe(paste("bedtools intersect -a ", file.score,  " -b ", file.refer, " -wa ", sep=" ")) ) );
	unlink(file.refer);
	unlink(file.score);
	return(ret);
}

bed_subtract <- function( score_bed, refer_bed ){
	file.refer <- tempfile(fileext=".bed");
	file.score <- tempfile(fileext=".bed");
	write.table( refer_bed[,c(1:3)], file=file.refer, quote=F, col.names=F, row.names=F, sep="\t");
	write.table( score_bed, file=file.score, quote=F, col.names=F, row.names=F, sep="\t");
	ret <- unique( read.table( pipe(paste("bedtools intersect -a ", file.score,  " -b ", file.refer, " -v ", sep=" ")) ) );
	unlink(file.refer);
	unlink(file.score);
	return(ret);
}


get.PWM<-function(seq.str)
{
   m <- matrix(0.01, ncol=4, nrow=nchar(seq.str[1]));
   for (i in 1:NROW(seq.str) )
   {
        seq.vec <- strsplit(seq.str[i], split="")[[1]];
        
        for( j in 1:NROW(seq.vec))
        {
            if (seq.vec[j]=="A") m[j,1]<-m[j,1]+1;
            if (seq.vec[j]=="C") m[j,2]<-m[j,2]+1;
            if (seq.vec[j]=="G") m[j,3]<-m[j,3]+1;
            if (seq.vec[j]=="T") m[j,4]<-m[j,4]+1;
        }
   }
   
   m <- m/rowSums(m);
   colnames(m)<- c("A", "C", "G", "T")
   
   return(data.frame(Pos=1:NROW(m), m) );
}

file.K562.DHS.peak  <- "./k562.merged.broad.peak.bed";
file.grocap="../k562/hg19.k562.new_hmm2b.post2.bed"

make_big_bed<-function( file.nonneg,dist=100 )
{
	system( paste("cat ../k562/wgEncodeOpenChromDnaseK562PkV2.narrowPeak ../k562/GSM646567_hg19_wgEncodeUwDgfK562Pk.narrowPeak.txt ../k562/GSM646567_hg19_wgEncodeUwDgfK562Pk.macs2.narrowPeak | awk -v OFS='\\t' '{print $1,$2,$3}' - > ", file.nonneg ), intern=TRUE );

	nonneg_bed <- unique(read.table( pipe(paste("cat ", file.nonneg, file.grocap, " | sort-bed ", file.nonneg," | bedtools merge -i - -d 100 ", sep=" " ) ) )[,c(1:3)]);
	nonneg_bed[,2] <- nonneg_bed[,2] - dist
	idx.mis <- which(nonneg_bed[,2]<0);
	if(length(idx.mis)>0) nonneg_bed[idx.mis,2] <- 0;
	nonneg_bed[,3] <- nonneg_bed[,3] + dist
	write.table( nonneg_bed, file=file.nonneg, quote=F, row.name=F, col.names=F, sep="\t" );
}

make_big_bed( file.K562.DHS.peak );


TATA.pwm <- get.PWM(TATA.seq);
PSE.pwm <- get.PWM(PSE.seq);

write.table(TATA.pwm, file="TATA.pwm", quote=F, col.names=T, row.names=F, sep="\t");
write.table(PSE.pwm, file="PSE.pwm", quote=F, col.names=T, row.names=F, sep="\t");

require(rtfbsdb);
db <- CisBP.extdata("Homo_sapiens");
tfs <- tfbs.createFromCisBP(db, motif_id="M5917_1.02");
tfs;

tfs <- tfbs.importMotifs(tfs, 'pwm.matrix', c( "TATA.pwm","PSE.pwm" ), c("TATA", "PSE"), header=TRUE );
show(tfs);

file.twoBit <- "/fs/cbsudanko/storage/data/hg19/hg19.2bit"
TRE.bed <- read.table("/local/workdir/zw355/proj/prj10-dreg/new-rf-201803/G1/G1.dREG.peak.score.bed.gz");

t1 <- tfbs.scanTFsite( tfs, file.twoBit,  TRE.bed, threshold = 1, ncores = 3);

tfbs.drawLogo(tfs, file.pdf="fig-S11-motif.pdf");


POLR2A.bed <- read.table(pipe("zcat /fs/cbsudanko/storage/data/hg19/all/ENCODE_tf_peak_calls/wgEncodeRegTfbsClusteredWithCellsV3.bed.gz | grep POLR2A - | grep  K562 "))
POLR2A.bed <- POLR2A.bed[,c(1:3)]
RPC155.bed <- read.table(pipe("zcat /fs/cbsudanko/storage/data/hg19/all/ENCODE_tf_peak_calls/wgEncodeRegTfbsClusteredWithCellsV3.bed.gz | grep RPC155 - | grep  K562 "))
RPC155.bed <- RPC155.bed[,c(1:3)]

DHS.bed <- read.table(file.K562.DHS.peak)
NONDHS.bed <- bed_subtract(TRE.bed, DHS.bed);
dREG.DHSplus.bed <- bed_intersect(TRE.bed, DHS.bed);
dREG.DHSminus.bed <- NONDHS.bed;

stat <- c();
for(i in 1:10)
{
     TATA.bed <- t1$result[[2]][ t1$result[[2]]$score>=i,];
     PSE.bed <- t1$result[[3]][ t1$result[[3]]$score>=i,];
     
     file.TATA <- tempfile(".bed");
     file.PSE <- tempfile(".bed");
     write.table( TATA.bed, file=file.TATA, quote=F, row.names=F, col.names=F, sep="\t");
     write.table( PSE.bed, file=file.PSE, quote=F, row.names=F, col.names=F, sep="\t");
     
     comb.bed <- read.table(pipe(paste("bedtools closest -a ", file.PSE, " -b ", file.TATA)))
     bed.Pol3 <- comb.bed[ which(comb.bed$V10-comb.bed$V3 < 40 & comb.bed$V10-comb.bed$V3 > 10 ),]
     bed.Pol3 <- bed.Pol3[,c(1,2,11)];
     
     bed.int.Pol2 <- bed.int.Pol3 <- bed.int.DHS <-c();
     if(NROW(bed.Pol3)>0)
     {
		 bed.int.Pol2 <- try(bed_intersect ( bed.Pol3, POLR2A.bed ))
		 bed.int.Pol3 <- try(bed_intersect ( bed.Pol3, RPC155.bed ))
		 bed.int.NONDHS <- try(bed_intersect ( bed.Pol3, NONDHS.bed ))
     }
     else
     	break;
    
     Pol3.fisher <- fisher.test(data.frame(c(NROW(RPC155.bed) - NROW(bed.int.Pol3), NROW(bed.int.Pol3)), c(NROW(NONDHS.bed)-NROW(bed.int.NONDHS), NROW(bed.int.NONDHS) ) ));
     Pol2.fisher <- fisher.test(data.frame(c(NROW(POLR2A.bed) - NROW(bed.int.Pol2), NROW(bed.int.Pol2)), c(NROW(NONDHS.bed)-NROW(bed.int.NONDHS), NROW(bed.int.NONDHS) ) ));
     Pol23.fisher <- fisher.test(data.frame(c(NROW(POLR2A.bed) - NROW(bed.int.Pol2), NROW(bed.int.Pol2)), c(NROW(RPC155.bed) - NROW(bed.int.Pol3), NROW(bed.int.Pol3) ) ));
     
     stat <- rbind(stat,  c( i, NROW(bed.int.Pol2), NROW(bed.int.Pol3), NROW(bed.int.NONDHS), NROW(POLR2A.bed), NROW(RPC155.bed), NROW(NONDHS.bed), p.Pol3=Pol3.fisher$p.value, p.Pol2=Pol2.fisher$p.value  ) )
}

colnames(stat) <- c( "score", "int.Pol2", "int.Pol3", "int.NONDHS", "bed.POLR2A", "bed.RPC155", "bed.DHS-", "pv.Pol3", "pv.Pol2" )
show(stat);
browser();

stat[,2] <- stat[,2] /NROW(POLR2A.bed);
stat[,3] <- stat[,3] /NROW(RPC155.bed);
stat[,4] <- stat[,4] /NROW(NONDHS.bed);


pdf("fig-S14.pdf")
cols <- c("#86db63", "#f92367", "#5070fb")
plot(1,1, type="n", xlim=c(1,10), ylim=c(0,0.055), xlab="Motif score", ylab="Fraction of peak")
points(1:9, stat[,2], pch=17, col=cols[1], cex=1)
lines (1:9, stat[,2], col=cols[1]);
points(1:9, stat[,3], pch=18, col=cols[2], cex=1)
lines (1:9, stat[,3], col=cols[2]);
points(1:9, stat[,4], pch=19, col=cols[3], cex=1)
lines (1:9, stat[,4], col=cols[3]);
legend("topright", c("Pol II w/ dREG", "Pol III w/ dREG", "DHS-"), text.col=cols, col=cols, pch=17:19, cex=1, ncol=1, horiz=F )
dev.off();

pdf("fig-S14-b.pdf")
barplot(cbind(PolII=NROW(POLR2A.bed), PolIII=NROW(RPC155.bed), DHS=NROW(NONDHS.bed)), col="gray");
dev.off();

