#$ -S /bin/bash
#$ -cwd
#$ -N tss_JURKAT
#$ -o tss_JURKAT.out.$JOB_ID
#$ -j y
#$ -pe bscb 16
#$ -M dankoc@gmail.com
#$ -m be
#$ -l h_rt=36:00:00
#$ -q long_term.q

STARTDIR=`pwd`

## Copy files to scratch space (/workdir and /SSD).
SCRATCH=/SSD/cgd24_tssDetector_cd4_Jurkat
mkdir $SCRATCH
cp /home/cgd24/projects/tss_detector/train_svm/asvm.intersDNase.getTrainSet.RData $SCRATCH ## 
cp /bscb/bscb07/cgd24/projects/CD4/Alignments/J-U-cc_*.bw $SCRATCH ## Combined
cp /bscb/bscb07/cgd24/projects/CD4/Alignments/J-PI-cc_*.bw $SCRATCH
cp /home/cgd24/projects/tss_detector/run_svm_cd4/scan_cd4-h1.R $SCRATCH
cd $SCRATCH

## Run R.
R --no-save < scan_cd4-h1.R
gzip *.predictions.bedGraph

## Copy data files back.
cp *.predictions.bedGraph.gz $STARTDIR
rm -Rf $SCRATCH
