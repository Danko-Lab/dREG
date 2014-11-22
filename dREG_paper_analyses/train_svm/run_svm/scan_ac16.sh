#$ -S /bin/bash
#$ -cwd
#$ -N tss_detector
#$ -o tss_detector.out.$JOB_ID
#$ -j y
#$ -pe bscb 32
#$ -M dankoc@gmail.com
#$ -m be
#$ -l h_rt=24:00:00

STARTDIR=`pwd`

## Copy files to scratch space (/workdir and /SSD).
SCRATCH=/SSD/cgd24_dREG_AC16
mkdir $SCRATCH

cp /home/cgd24/projects/tss_detector/run_svm_ac16/scan_ac16.R $SCRATCH ## 
cp /home/cgd24/projects/tss_detector/train_svm/asvm.intersDNase.getTrainSet.RData $SCRATCH ## 
cp /bscb/bscb07/cgd24/data/hg19/ac16/groseq/ac16*.bw $SCRATCH ## bigWig files.

cd $SCRATCH

## Run R.
R --no-save < scan_ac16.R
gzip ac16.predictions.bedGraph

## Copy data files back.
cp ac16.predictions.bedGraph.gz $STARTDIR
cp ac16.predictions.bedGraph.gz ~/nextgen/home/cgd24/work/tss_detector/

rm -rf $SCRATCH
