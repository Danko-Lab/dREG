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
SCRATCH=/SSD/cgd24_tssDetector_cd4
mkdir $SCRATCH
cp /home/cgd24/projects/tss_detector/run_svm_cd4/scan_cd4-h1.R $SCRATCH ## 
cp /home/cgd24/projects/tss_detector/train_svm/asvm.RData $SCRATCH ## 
cp ~/nextgen/projects/GROseq/CD4/Alignments/H1-U_*.bw $SCRATCH
cp ~/nextgen/projects/GROseq/CD4/Alignments/J-U_*.bw $SCRATCH 
cp ~/nextgen/projects/GROseq/CD4/Alignments/J-U-cc_*.bw $SCRATCH
cp ~/nextgen/projects/GROseq/CD4/Alignments/H1-PI_*.bw $SCRATCH
cp ~/nextgen/projects/GROseq/CD4/Alignments/J-PI_*.bw $SCRATCH 
cp ~/nextgen/projects/GROseq/CD4/Alignments/J-PI-cc_*.bw $SCRATCH
cd $SCRATCH

## Run R.
R --no-save < scan_cd4-h1.R
gzip *.predictions.bedGraph

## Copy data files back.
cp *.predictions.bedGraph.gz $STARTDIR
rm -Rf $SCRATCH
