#$ -S /bin/bash
#$ -cwd
#$ -N tss_detector
#$ -o tss_detector.out.$JOB_ID
#$ -j y
#$ -pe bscb 16
#$ -M dankoc@gmail.com
#$ -m be
#$ -l h_rt=48:00:00
#$ -q long_term.q

STARTDIR=`pwd`

## Copy files to scratch space (/workdir and /SSD).
SCRATCH=/SSD/cgd24_dREG_mcf7_40
mkdir $SCRATCH

cp /home/cgd24/projects/tss_detector/run_svm_MCF7/scan_mcf7.40m.R $SCRATCH ## 
cp /home/cgd24/projects/tss_detector/train_svm/asvm.intersDNase.getTrainSet.RData $SCRATCH ## 
cp /bscb/bscb07/cgd24/data/hg19/mcf7/groseq/*.bw $SCRATCH ## bigWig files.

cd $SCRATCH

## Run R.
R --no-save < scan_mcf7.40m.R
gzip mcf7.E2_40m.predictions.bedGraph

## Copy data files back.
cp mcf7.E2_40m.predictions.bedGraph.gz $STARTDIR

rm -rf $SCRATCH
