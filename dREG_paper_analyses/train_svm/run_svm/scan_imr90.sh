#$ -S /bin/bash
#$ -cwd
#$ -N tss_detector
#$ -o tss_detector.out.$JOB_ID
#$ -j y
#$ -pe bscb 16
#$ -M dankoc@gmail.com
#$ -m be
#$ -l h_rt=24:00:00

STARTDIR=`pwd`

## Copy files to scratch space (/workdir and /SSD).
SCRATCH=/SSD/cgd24_dREG_IMR90
mkdir $SCRATCH

cp /home/cgd24/projects/tss_detector/run_svm_IMR90/scan_imr90.R $SCRATCH ## 
cp /home/cgd24/projects/tss_detector/train_svm/asvm.intersDNase.getTrainSet.RData $SCRATCH ## 
cp /bscb/bscb07/cgd24/data/hg19/imr90/groseq/groseq_*.bigWig $SCRATCH ## bigWig files.

cd $SCRATCH

## Run R.
R --no-save < scan_imr90.R
gzip imr90.predictions.bedGraph

## Copy data files back.
cp imr90.predictions.bedGraph.gz $STARTDIR

rm -rf $SCRATCH
