#$ -S /bin/bash
#$ -cwd
#$ -N dREG.hela
#$ -o dREG.hela.$JOB_ID
#$ -j y
#$ -pe bscb 32
#$ -M dankoc@gmail.com
#$ -m be
#$ -l h_rt=24:00:00

STARTDIR=`pwd`

## Copy files to scratch space (/workdir and /SSD).
SCRATCH=/SSD/cgd24_dREG_HELA
mkdir $SCRATCH

cp /home/cgd24/projects/tss_detector/run_svm_hela/scan_hela.R $SCRATCH ## 
cp /home/cgd24/projects/tss_detector/train_svm/asvm.RData $SCRATCH ## 
cp ~/nextgen/data/GROseq.parser/hg19/hela/groseq/HeLa*.bw $SCRATCH ## bigWig files.

cd $SCRATCH

## Run R.
R --no-save < scan_hela.R
gzip hela.predictions.bedGraph

## Copy data files back.
cp hela.predictions.bedGraph.gz $STARTDIR

rm -rf $SCRATCH
