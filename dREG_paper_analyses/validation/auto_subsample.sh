#!/bin/bash

THRESHOLD=0.7
CHINFO=/usr/data/GROseq.parser/hg19/k562/proseq/chromInfo.hg19

# The input are subsample numbers, in units of 1 million.

# Checking to see if we've already subsampled.
#to_subsample=()
#for n in "$@"
#do
#    if [ ! -e /usr/projects/GROseq.parser/tss_detecter/subsample/K562_unt.sort.subsamp_"$n"m.bed ]
#        then
#            echo "Missing .bed file for "$n
#            to_subsample +=($n)
#    fi

# Creating necessary subsampled files.
echo "Subsampling."
python subsample.py $@
#python subsample.py $to_subsample

# Convert them to bigwig. (this is just a modification of Charles' bedToBigWig)
echo "Converting to bigWig."
for n in "$@"
do
    #gzip /usr/projects/GROseq.parser/tss_detecter/subsample/K562_unt.sort.subsamp_"$n"m.bed
    THEFILE=/usr/projects/GROseq.parser/tss_detecter/subsample/K562_unt.sort.subsamp_"$n"m.bed
    echo $THEFILE

    ## Remove rRNA and reverse the strand (PRO-seq).
    cat $THEFILE | grep "rRNA" -v | grep "_" -v | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,$6=="+"?"-":"+"}' | gzip > $THEFILE.nr.rs.bed.gz
    
    ## Convert to wig
    /usr/data/GROseq.parser/hg19/k562/proseq/bedToWig $THEFILE.nr.rs.bed.gz $THEFILE
    gzip $THEFILE*.wig

    ## Then to bigWig
    wigToBigWig $THEFILE\_plus.wig.gz $CHINFO $THEFILE.gz_plus.bw
    wigToBigWig $THEFILE\_minus.wig.gz $CHINFO $THEFILE.gz_minus.bw

    rm $THEFILE.nr.rs.bed.gz
    echo "Running test_train_svm_subsampled.R..."
    R --slave --vanilla --file=/home/sh985/featureDetector/test_functions/test_train_svm_subsampled.R --args $n\m $THRESHOLD
done
