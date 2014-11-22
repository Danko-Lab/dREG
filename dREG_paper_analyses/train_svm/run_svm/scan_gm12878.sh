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
SCRATCH=/SSD/cgd24_dREG_gm12878
mkdir $SCRATCH
cp /home/cgd24/projects/tss_detector/train_svm/scan_gm12878.R $SCRATCH ## 
#cp /home/cgd24/projects/tss_detector/train_svm/asvm.RData $SCRATCH ## 
cp /home/cgd24/projects/tss_detector/train_svm/asvm.intersDNase.getTrainSet.RData $SCRATCH
cp /home/cgd24/projects/tss_detector/data/gm12878/groseq_*.bigWig $SCRATCH ## bigWig files.
cd $SCRATCH

## Run R.
R --no-save < scan_gm12878.R
gzip gm12878.predictions.bedGraph

## Copy data files back.
cp gm12878.predictions.bedGraph.gz $STARTDIR

rm -Rf $SCRATCH
