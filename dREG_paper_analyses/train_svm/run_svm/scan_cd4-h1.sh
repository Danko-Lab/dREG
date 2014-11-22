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
SCRATCH=/SSD/cgd24_tssDetector_cd4_Jurkat
mkdir $SCRATCH
cp /home/cgd24/projects/tss_detector/train_svm/asvm.intersDNase.getTrainSet.RData $SCRATCH ## 
cp /bscb/bscb07/cgd24/projects/CD4/Alignments/J-U-cc_*.bw $SCRATCH ## Combined
cp /bscb/bscb07/cgd24/projects/CD4/Alignments/J-PI-cc_*.bw $SCRATCH
cd $SCRATCH

## Run R.
R --no-save < scan_cd4-h1.R
gzip *.predictions.bedGraph

## Copy data files back.
cp *.predictions.bedGraph.gz $STARTDIR
rm -Rf $SCRATCH
