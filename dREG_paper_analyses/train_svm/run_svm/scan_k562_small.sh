#$ -S /bin/bash
#$ -cwd
#$ -N dREG.k562_small
#$ -o dREG.k562_small.$JOB_ID
#$ -j y
#$ -pe bscb 30
#$ -M dankoc@gmail.com
#$ -m be
#$ -l h_rt=24:00:00

STARTDIR=`pwd`

## Copy files to scratch space (/workdir and /SSD).
SCRATCH=/SSD/cgd24_dREG_k562_small
mkdir $SCRATCH

cp /home/cgd24/projects/tss_detector/run_svm_k562_small/scan_k562_small.R $SCRATCH ## 
cp /home/cgd24/projects/tss_detector/train_svm/asvm.RData $SCRATCH ## 
cp ~/nextgen/data/GROseq.parser/hg19/k562/groseq/groseq_*.bigWig $SCRATCH ## bigWig files.

cd $SCRATCH

## Run R.
R --no-save < scan_k562_small.R
gzip k562_small.predictions.bedGraph

## Copy data files back.
cp k562_small.predictions.bedGraph.gz $STARTDIR

rm -rf $SCRATCH
