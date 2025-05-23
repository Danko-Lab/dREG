dREG
===============

Detection of Regulatory DNA Sequences using GRO-seq Data.

Cloud Computing Service:
-------------------------
We provide a cloud service that allows users to run dREG on a GPU server:

https://dreg.dnasequence.org/

![Hi](https://github.com/Danko-Lab/dREG/raw/master/dreg-gateway.png?v=4&s=200 "dREG gateway")

### Server status

The dREG gateway runs on the [Bridges2 supercomputer](https://www.psc.edu/resources/bridges-2/) at the Pittsburgh Supercomputing Center, which is in turn supported by the [NSF ACCESS](https://access-ci.org/) program. If you're having issues with the gateway, please check the [ACCESS server status updates log](https://operations.access-ci.org/infrastructure_news) for any ongoing Bridges2 service outages.

### Important note for the Exchange email users:

The Exchange email system might quarantine all emails including the word  “password” or other sensitive stuffs in links. (https://technet.microsoft.com/en-us/library/aa997692(v=exchg.160).aspx).

Unfortunately, some emails from dREG gateway are quarantined by this spam policy. Usually, these quarantined emails are not delivered to the email box, so they cannot be checked in any email folders, including junk, spam or inbox. If you find that emails from dREG gateway are not delivered into your email box, please conect the administrator of your email system. For the Cornell email, please check this link:

https://it.cornell.edu/spam-control/log-quarantine-management-spam-control

Abstract:
--------
Identification of the genomic regions that regulate transcription remains an important open problem.  We have recently shown that global run-on and sequencing (GRO-seq) with enrichment for 5-prime-capped RNAs reveals patterns of divergent transcription that accurately mark active transcriptional regulatory elements (TREs), including enhancers and promoters.  Here, we demonstrate that active TREs can be identified with comparable accuracy by applying sensitive machine-learning methods to standard GRO-seq and PRO-seq data, allowing TREs to be assayed together with transcription levels, elongation rates, and other transcriptional features, in a single experiment.  Our method, called discriminative Regulatory Element detection from GRO-seq (dREG), summarizes GRO-seq read counts at multiple scales and uses support vector regression to predict active TREs.  The predicted TREs are strongly enriched for marks associated with functional elements, including H3K27ac, transcription factor binding sites, eQTLs, and GWAS-associated SNPs.  Using dREG, we survey TREs in eight cell types and provide new insights into global patterns of TRE assembly and function. 

Data preparation: 
==========================

dREG takes bigWig files with double strands as the input. The bigWig files should follow 3 rules:

1) Each read is mapped at 5’ (GRO-seq) or 3’ (PRO-seq) position (point mode) , not mapped to a continuous region starting from 5’ or 3’.  This is different with the software Tfit.

2) Only positive values or only negative values in each strand, no mixture.

3) No normalization

As for how to generate bigWig files from fastq data, please refer to https://github.com/Danko-Lab/proseq2.0/.


Installation instructions: 
==========================

dREG can either be installed from source as described below, or via a [docker image of dREG](https://hub.docker.com/r/biohpc/dreg) created by the BioHPC staff at Cornell University.

Supported OS:
-------------
Linux and Mac OSX are currently supported.

Required software
-----------------
* bedops (http://bedops.readthedocs.org/en/latest/index.html)
* R (http://www.r-project.org/)
* bigWig R package (https://github.com/andrelmartins/bigWig).
* tabix (frmo htslib)
* bedgraphtobigwig (from ucsc tools)

This software is already installed on many UNIX systems.  Users can install the most appropriate version of these files for Ubuntu using: 

    sudo apt-get install r-base-core
    sudo apt-get install libssl1.0.0 libssl-dev

Users who are not sure how to install the proper dependencies on their system should check with their system administrator for help.  

dREG also has several dependencies within R.  These include **data.table**, **e1071**, **mvtnorm**, **parallel**, **randomForest**, **rmutil**, **rphast**, and **snowfall**.  These packages are all availiable on the CRAN repository.  For convenience, users can install these packages using the makefile:

    make R_dependencies

If users run into any problems, they should contact the package author for assistance.

Install dREG
------------
Users should change to the directory containing this README.md file, and can then install dREG by typing the following:

(1a) Install R dependencies

    make R_dependencies

(1b) Install dREG

    make dreg

Download dREG models
--------------------
The most recent pre-trained dREG model can be downloaded from one of the following sources:

(1) https://dreg.dnasequence.org/themes/dreg/assets/file/asvm.gdm.6.6M.20170828.rdata

(2) ftp://cbsuftp.tc.cornell.edu/danko/hub/dreg.models/asvm.gdm.6.6M.20170828.rdata

(3) https://zenodo.org/records/10113379

If you have issues downloading this model file, please contact us.

Usage instructions:
===================

dREG provides two solutions to identify TREs in this R package. 

The first solution implemented in the early package, is to ***predict dREG scores*** and detect the broad dREG peaks with the aid of Perl program. In order to identify narrow peak, these broad peaks need to be refined using [dREG-HD package](https://github.com/Danko-Lab/dREG.HD). This solution is generally no longer needed, but we have maintained this functionality in the package for backwards compatibility with older dREG SVR models.

The second solution implements the ***peak calling*** function using the dREG scores based on the improved SVR model. Compared with the broad peaks in the first solution, this solution directly generates the narrow peaks with peak score, probability, center position. **This is the primary method that we recommend, as it provides an all-in-one peak calling method**. However, this is only compatible with the newer dREG SVR model, which requires GPU acceleration. If you do not have access to GPUs, please try our online computational gateway described above. 

## 1) Peak calling (recommended)

To use this solution, type: 

    bash run_dREG.bsh plus_strand.bw minus_strand.bw out_prefix asvm.RData [nthreads] [GPU]

    plus_strand.bw	-- PRO-seq data (plus strand).  Read counts (not normalized) formatted as a bigWig file.
    minus_strand.bw	-- PRO-seq data (minus strand). Read counts (not normalized) formatted as a bigWig file.
    out_prefix      -- The prefix of the output file.
    asvm.RData      -- The path to the RData file containing the pre-trained SVM.
    cpu_cores       -- [optional, default=1] indicating how many CPU cores can be used.
    gpu_id          -- [optional, default=NA] indicating GPU id when multiple GPU cards are available.
    

For example, to run dREG on the PRO-seq data, use:

    bash run_dREG.bsh proseq.plus.bw proseq.minus.bw proseq.test asvm.gdm.6.6M.20170828.rdata 16 0

Three files below are generated in this solution:  

1. <out_prefix>.dREG.infp.bed.gz:

    all informative sites with dREG scores.

2. <out_prefix>.dREG.peak.full.bed.gz: 

    full peak information, including peak position, max score, probability, center.

3. <out_prefix>.dREG.peak.score.bed.gz:

    reduced peak information, only including peak position, max score.

**Notice:** 

(1) This solution doesn't work with the model trained before 2017. The new SVR model can be downloaded as described above.

(2) That command takes 4~12 hours to execute on NVIDA K80 GPU using Rgtsvm package. Due to very long computational time, we don't suggest to run peak calling on CPU, even in parallel mode.

## 2) Predicting dREG scores (legacy)

For this solution, dREG takes three files as input, and outputs one file.  Input files include the PRO-seq read distributions on the plus and minus strand (which are separate files), and parameters of the pre-trained support vector regression (SVR) model.  

* PRO-seq files are required to be in the bigWig format standard created by the UCSC (more information can be found here: http://genome.ucsc.edu/goldenPath/help/bigWig.html).  
* The SVR model is included in this package (under dREG_model/asvm.RData).  Users are advised to use that when possible.

To use dREG, type: 

    bash run_predict.bsh plus_strand.bw minus_strand.bw out_prefix asvm.RData [nthreads] [GPU]

    plus_strand.bw	-- PRO-seq data (plus strand).  Read counts (not normalized) formatted as a bigWig file.
    minus_strand.bw	-- PRO-seq data (minus strand). Read counts (not normalized) formatted as a bigWig file.
    out_prefix      -- The prefix of the output file.
    asvm.RData      -- The path to the RData file containing the pre-trained SVM.
    cpu_cores       -- [optional, default=1] indicating how many CPU cores can be used.
    gpu_id          -- [optional, default=NA] indicating GPU id when multiple GPU cards are available.

For example, to run dREG on the example data (PRO-seq from chr21 in K562 cells), use:

    bash run_predict.bsh example/K562.chr21.plus.bw example/K562.chr21.minus.bw k562.test dREG_model/asvm.RData 2 

That command takes ~2-3 hours to execute on Ubuntu on a core i5 desktop computer (CPU only).

If GPU is available with 16 CPU cores, use:

    bash run_predict.bsh example/K562.chr21.plus.bw example/K562.chr21.minus.bw k562.test dREG_model/asvm.RData 16 1

dREG outputs a bedGraph file of scores.  If desired, users can convert this file into a merged file of dREG 'peaks', or regions which fit the profile of a transcribed regulatory element.   For convenience, users can use the included bash script (writeBed.bsh) to identify dREG peaks.  This script is used as follows:

    bash writeBed.bsh 0.8 out_prefix.bedGraph.gz

Here `0.8` denotes the threshold to call a regulatory element, and the out_prefix.bedGraph.gz is the output of the dREG run.  Note that this feature requires the bedOps package (https://bedops.readthedocs.org/en/latest/).

The threshold `0.8` is used to the predictions from the SVR models in this page(https://github.com/Danko-Lab/dREG-Model). For the huge SVR models, we suggest to use `0.25` as threshold in this solution. 


# Document

dREG is an R package, and that provides some additional flexibility for users familiar with R. Currently you can get details about each R function from the dREG manual (https://github.com/Danko-Lab/dREG/blob/master/dREG-manual.pdf).  We are actively working to document each function in the package.  

How to cite
===================
(1) Danko, C. G., Hyland, S. L., Core, L. J., Martins, A. L., Waters, C. T., Lee, H. W., ... & Siepel, A. (2015). Identification of active transcriptional regulatory elements from GRO-seq data. Nature methods, 12(5), 433-438. (https://www.nature.com/articles/nmeth.3329)

(2) Wang, Z., Chu, T., Choate, L. A., & Danko, C. G. (2018). Identification of regulatory elements from nascent transcription using dREG. bioRxiv, 321539. (https://www.biorxiv.org/content/early/2018/05/14/321539.abstract)

