#
# Sets up a run of PIQ on the same data as the naive method.

require(rtfbsdb)


NRF1 <- read.table("/usr/data/GROseq.parser/pwm_data/jolma/teal/NRF1.YGCGCATGCGCN.pwm", header=TRUE)
ELF1 <- read.table("/usr/data/GROseq.parser/pwm_data/jolma/teal/ELF1.AACCCGGAAGTR.pwm", header=TRUE)
SP1  <- read.table("/usr/data/GROseq.parser/pwm_data/jolma/teal/SP1.GCCMCGCCCMC.pwm", header=TRUE)
MAX  <- read.table("/usr/data/GROseq.parser/pwm_data/jolma/teal/MAX.NNCACGTGNN.pwm", header=TRUE)

write.table(t(NRF1), "jaspar.format.txt")
write.table(t(ELF1), "jaspar.format.txt", append=TRUE)
write.table(t(SP1), "jaspar.format.txt", append=TRUE)
write.table(t(MAX), "jaspar.format.txt", append=TRUE)

BAM <- "/usr/data/GROseq.parser/hg19/k562/dnase/GSM646567_hg19_wgEncodeUwDgfK562Aln.bam"

