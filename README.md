dREG
===============

Detection of Regulatory DNA Sequences using GRO-seq Data.

Online computational gateway
-------------------------
We provide a computational gateway to run dREG on GPU server, the users don't need to install any software, only upload the bigWig files and wait for the results, it is simple and easy. Please click the link to try this site:

https://dreg.dnasequence.org/

### Important note for the Exchange email users:

The Exchange email system might quarantine all emails including the word  “password” or other sensitive stuffs in links. (https://technet.microsoft.com/en-us/library/aa997692(v=exchg.160).aspx).

Unfortunately, some emails from dREG gateway are quarantined by this spam policy. Usually these quarantined emails are not delivered to the email box, so they can not be checked in any email folders, including junk, spam or inbox. If you find the emails from dREG gateway are not delivered into your email box, please conect the administrator of your email system. For the Cornell email, please check this link:

https://it.cornell.edu/spam-control/log-quarantine-management-spam-control

Abstract
--------
Identification of the genomic regions that regulate transcription remains an important open problem.  We have recently shown that global run-on and sequencing (GRO-seq) with enrichment for 5-prime-capped RNAs reveals patterns of divergent transcription that accurately mark active transcriptional regulatory elements (TREs), including enhancers and promoters.  Here, we demonstrate that active TREs can be identified with comparable accuracy by applying sensitive machine-learning methods to standard GRO-seq and PRO-seq data, allowing TREs to be assayed together with transcription levels, elongation rates, and other transcriptional features, in a single experiment.  Our method, called discriminative Regulatory Element detection from GRO-seq (dREG), summarizes GRO-seq read counts at multiple scales and uses support vector regression to predict active TREs.  The predicted TREs are strongly enriched for marks associated with functional elements, including H3K27ac, transcription factor binding sites, eQTLs, and GWAS-associated SNPs.  Using dREG, we survey TREs in eight cell types and provide new insights into global patterns of TRE assembly and function. 


Installation instructions: 
==========================

dREG will ultimately be availiable in the R repository CRAN to ease installation, and source code will be availiable on GitHub (https://github.com/Danko-Lab/dREG).  

Supported OS:
-------------
Linux and Mac OSX are currently supported.

Required software
-----------------
* R (http://www.r-project.org/)
* mysql-dev (http://dev.mysql.com/downloads/).
* bigWig R package (https://github.com/andrelmartins/bigWig; will be public very soon).

This software is already installed on many UNIX systems.  Users can install the most appropriate version of these files for Ubuntu using: 

    sudo apt-get install r-base-core
    sudo apt-get install libmysqlclient-dev
    sudo apt-get install libssl1.0.0 libssl-dev

Users who are not sure how to install the proper dependencies on their system should check with their system administrator for help.  

dREG also has several dependencies within R.  These include rphast, grid, boot, e1071, and parallel.  These packages are all availiable on the CRAN repository.  For convenience, users can install these packages using the makefile:

    make R_dependencies

If users run into any problems they should contact the package author for assistance.

Install dREG
------------
Users should change to the directory containing this README.md file, and can then install dREG by typing the following:

(1a) Install R dependencies

    make R_dependencies

(1b) Install dREG

    make dreg

Get the dREG models
-------------------
Pre-trained models that can be used to predict dREG scores across the genome are availiable in mammals and drosophila.  Get the appropriate model for your system here: https://github.com/Danko-Lab/dREG-Model

or download the newest model from FTP:
ftp://cbsuftp.tc.cornell.edu/danko/hub/dreg.models/asvm.gdm.6.6M.20170828.rdata

Usage instructions:
===================

dREG provides two solutions to identify TREs in this R package. 

The first solution implemented in the early package, is to ***predict dREG scores*** and detect the broad dREG peaks with the aid of Perl program. In order to identify narrow peak, these broad peaks need to be refined using dREG-HD package.

The second solution implements the peak calling function using the dREG scores based on the imporved SVR model. Compared with the broad peaks in the first solution, this solution generates the narrow peaks with peak score, probability, center position.

## 1) Predicting dREG scores

For this solution, dREG takes three files as input, and outputs one file.  Input files include the PRO-seq read distributions on the plus and minus strand (which are separate files), and parameters of the pre-trained support vector regression (SVR) model.  

* PRO-seq files are required to be in the bigWig format standard created by the UCSC (more information can be found here: http://genome.ucsc.edu/goldenPath/help/bigWig.html).  
* The SVR model is included in this package (under dREG_model/asvm.RData).  Users are advised to use that when possible.

To use dREG, type: 

    bash run_dREG.bsh plus_strand.bw minus_strand.bw out_prefix asvm.RData [nthreads] [GPU]

    plus_strand.bw	-- PRO-seq data (plus strand).  Read counts (not normalized) formatted as a bigWig file.
    minus_strand.bw	-- PRO-seq data (minus strand). Read counts (not normalized) formatted as a bigWig file.
    out_prefix      -- The prefix of the output file.
    asvm.RData      -- The path to the RData file containing the pre-trained SVM.
    [nthreads]      -- [optional, default=1] The number of threads to use.
    [GPU]           -- [optional, GPU or _blank_, default=_blank_] GPU can be used in this operation through the Rgtsvm package.


For example, to run dREG on the example data (PRO-seq from chr21 in K562 cells), use:

    bash run_dREG.bsh example/K562.chr21.plus.bw example/K562.chr21.minus.bw k562.test dREG_model/asvm.RData 2

If GPU is available with 16 CPU cores, use:

    bash run_dREG.bsh example/K562.chr21.plus.bw example/K562.chr21.minus.bw k562.test dREG_model/asvm.RData 15 GPU

That command takes ~2-3 hours to execute on Ubuntu on a core i5 desktop computer (CPU version).

dREG outputs a bedGraph file of scores.  If desired, users can convert this file into a merged file of dREG 'peaks', or regions which fit the profile of a transcribed regulatory element.   For convenience, users can use the included bash script (writeBed.bsh) to identify dREG peaks.  This script is used as follows:

    bash writeBed.bsh 0.8 out_prefix.bedGraph.gz

Here `0.8` denotes the threshold to call a regulatory element, and the out_prefix.bedGraph.gz is the output of the dREG run.  Note that this feature requires the bedOps package (https://bedops.readthedocs.org/en/latest/).

The threshold `0.8` is used to the predictions from the SVR models in this page(https://bedops.readthedocs.org/en/latest/). For the huge SVR models, we suggest to use `0.25` as threshold in this solution. 

## 2) Peak calling

To use this solution, type: 

    bash run_peakcalling.bsh plus_strand.bw minus_strand.bw out_prefix asvm.RData [nthreads] [GPU]

    plus_strand.bw	-- PRO-seq data (plus strand).  Read counts (not normalized) formatted as a bigWig file.
    minus_strand.bw	-- PRO-seq data (minus strand). Read counts (not normalized) formatted as a bigWig file.
    out_prefix      -- The prefix of the output file.
    asvm.RData      -- The path to the RData file containing the pre-trained SVM.
    [nthreads]      -- [optional, default=1] The number of threads to use.
    [GPU]           -- [`suggested`, GPU or _blank_, default=_blank_] GPU can be used in this operation through the Rgtsvm package.

For example, to run dREG on the PRO-seq data, use:

    bash run_peakcalling.bsh proseq.plus.bw proseq.minus.bw proseq.test asvm.gdm.6.6M.20170828.rdata 16 GPU

Three files below are generated in this solution:  

1. <out_prefix>.dREG.infp.bed.gz
   all informative sites with dREG scores.

2. <out_prefix>.dREG.peak.full.bed.gz 
   full peak information, including peak position, max score, probability, center.

3. <out_prefix>.dREG.peak.score.bed.gz
   reduced peak information, only including peak position, max score.

That command takes 8~12 hours to execute on NVIDA K80 GPU using Rgtsvm package.

# Reference

dREG is an R package, and that provides some additional flexibility for users familiar with R. Currently you can get details about each R function from the dREG manual (https://github.com/Danko-Lab/dREG/blob/master/dREG-manual.pdf).  We are actively working to document each function in the package.  
