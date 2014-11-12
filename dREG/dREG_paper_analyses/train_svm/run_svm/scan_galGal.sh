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
SCRATCH=/SSD/cgd24_tssDetector_galGal
mkdir $SCRATCH
cp /home/cgd24/projects/tss_detector/run_svm_galGal/scan_galGal.R $SCRATCH ## 
cp /home/cgd24/projects/tss_detector/train_svm/asvm.RData $SCRATCH ## 
cp ~/nextgen/home/cgd24/work/tss_detector/galGal/*.bw $SCRATCH
cd $SCRATCH

## Run R.
R --no-save < scan_galGal.R
gzip *.bedGraph

## Copy data files back.
cp *.bedGraph.gz $STARTDIR
cp *.bedGraph.gz ~/nextgen/home/cgd24/work/tss_detector/galGal/
rm -Rf $SCRATCH
