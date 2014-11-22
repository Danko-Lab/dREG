#$ -S /bin/bash
#$ -cwd
#$ -N tss_detector
#$ -o tss_detector.out.$JOB_ID
#$ -j y
#$ -pe bscb 16
#$ -M dankoc@gmail.com
#$ -m be
#$ -l h_rt=96:00:00
#$ -q long_term.q

STARTDIR=`pwd`

## Copy files to scratch space (/workdir and /SSD).
SCRATCH=/SSD/cgd24_scan_k562/
mkdir $SCRATCH
cp /home/cgd24/projects/tss_detector/train_svm/scan_k562.R $SCRATCH ## 
#cp /home/cgd24/projects/tss_detector/train_svm/asvm.RData $SCRATCH ## 
cp /home/cgd24/projects/tss_detector/train_svm/asvm.intersDNase.getTrainSet.RData $SCRATCH
cp /home/cgd24/projects/tss_detector/data/k562/K562_unt.sort.bed.gz_*.bw $SCRATCH ## bigWig files.
cd $SCRATCH

## Run R.
R --no-save < scan_k562.R
gzip k562.predictions.dn.bedGraph

## Copy data files back.
cp k562.predictions.dn.bedGraph.gz $STARTDIR

rm -rf $SCRATCH
