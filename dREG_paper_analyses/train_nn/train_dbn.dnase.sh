#$ -S /bin/bash
#$ -cwd
#$ -N tss_detector
#$ -o tss_detector.out.$JOB_ID
#$ -j y
#$ -pe bscb 6
#$ -M dankoc@gmail.com
#$ -m be
#$ -l h_rt=72:00:00

STARTDIR=`pwd`

## Copy files to scratch space (/workdir and /SSD).
SCRATCH=/SSD/cgd24_traindbn_dnase/
mkdir $SCRATCH
cd $SCRATCH

cp /home/cgd24/projects/tss_detector/train_dbn/train_dbn.dnase.R $SCRATCH ## 

## K562
cp /home/cgd24/projects/tss_detector/andre_hmm/hg19.k562.new_hmm2b.post2.bed $SCRATCH ## Andre's HMM Predictions
cp /home/cgd24/projects/tss_detector/data/GencodeMerge.IntersectOpStrand.bed $SCRATCH ## Gene overlap files
cp /home/cgd24/projects/tss_detector/data/k562/chromHmm.k562.enh.prom.bed.gz $SCRATCH ## Ernst chromHMM.
cp /home/cgd24/projects/tss_detector/data/k562/K562_unt.sort.bed.gz_*.bw $SCRATCH ## bigWig files.
cp /home/cgd24/nextgen/data/GROseq.parser/hg19/k562/dnase/wgEncodeUWDukeDnaseK562.fdr01peaks.hg19.bed $SCRATCH ## DNAse-1

## CD4
cp /home/cgd24/nextgen/data/GROseq.parser/hg19/cd4/dnase1fp/dnase1.peaks_peaks.narrowPeak $SCRATCH ## DNAse-1
cp ~/nextgen/projects/GROseq/NHP/AllData/All_Merge/H*.bw $SCRATCH ## PRO-seq
cp /home/cgd24/nextgen/data/GROseq.parser/hg19/cd4/chromhmm/CD4.chromHMM.Ernst2010.hg19.Prom.Enh.bed $SCRATCH ## Chromatin marks.

## HELA
cp /home/cgd24/nextgen/data/GROseq.parser/hg19/hela/groseq/HeLa*.bw $SCRATCH ## GRO-seq
cp /home/cgd24/nextgen/data/GROseq.parser/hg19/hela/dnase/uw.merge.narrowPeak.bed $SCRATCH ## DNAse-1
cp /home/cgd24/nextgen/data/GROseq.parser/hg19/hela/chromhmm/helas3.chromhmm.prom.enh.ins.bed $SCRATCH ## chromHMM

## MCF-7
cp /home/cgd24/nextgen/data/GROseq.parser/hg19/mcf7/groseq/MCF7.unt.all_*.bw $SCRATCH ## GRO-seq
cp /home/cgd24/nextgen/data/GROseq.parser/hg19/mcf7/dnase/wgEncodeAwgDnaseUwdukeMcf7UniPk.narrowPeak.gz $SCRATCH ## DNAse-1
cp /home/cgd24/nextgen/data/GROseq.parser/hg19/mcf7/histones/H3K4me1.liftOver.peaks.bed.gz $SCRATCH ## H3K4me1 (filling in for Ernst classes).

## Run R.
R --no-save < train_dbn.dnase.R

## Copy data files back.
cp *.bedGraph $STARTDIR
cp *.RData $STARTDIR

## Cleanup
rm -Rf $SCRATCH
