## ---------------------------
##
## Script name: 20220420_A_B_C_job-rscript10.R
##
## Purpose: This is an R script used to execute DADA2 denoising alogrithim on demultiplexed paired-end Illumina sequencing data. See pre-print: Paulson A.R., Lougheed, S.C., Huang, D., and Colautti, R.I. 2022. Evidence for symbionts and pathogens of blacklegged ticks (<i>Ixodes scapularis</i>) in an emerging Lyme disease hotspot. ##Preprint - webpage link to be added##
##
## Author: Amber R Paulson, Ph.D.
##
## Date Created: 2021-12-06
##
##
## ---------------------------
##
## Notes: 
## 
## See the READ.ME for more information on how to get started - https://github.com/damselflywingz/tick_microbiome
##
## Sep 5, 21 updates: remove primers, me2 (default) and triming for quality, dada2 = pooled for inference, added collapse stage back, repaired the tracker, increased max converge to 20, changed right to 239 trimm left, also truncation of the merged sequqnce from 250-265 bp
## April 20, 2022: added script details, section dividers and removed unecessary lines  
##     
##
## ---------------------------



##-------------------------------------------------------------------------
##                    Step 1 - R packages                                 -
##-------------------------------------------------------------------------

## below can be used to run R module on Queen's Universities Centre for Advanced Computation high-performance computing clusters but do not use if you are submitting this script to the SLURM scheduler

#job-rscript.R 
#module load StdEnv/2020 r/4.1.0
#R

## to exit the R module
#q() to quit

## all packages must be installed before submitting the R script to SLURM, see below for installing packages on the HPC cluster.

#install.packages("BiocManager", repos="https://mirror.rcg.sfu.ca/mirror/CRAN/")
#BiocManager::install(version = "3.13")

## update to the correct pathway where the packages need to be installed
#install.packages("/global/home/hpc4605/R/x86_64-pc-linux-gnu-library/4.1", repos = NULL, type = "win.binary")
#library(RCurl)

#BiocManager::install("dada2")
library('dada2')
#install.packages("here", repos="https://mirror.rcg.sfu.ca/mirror/CRAN/")
library('here')

##-------------------------------------------------------------------------
##                    Step 2 - Filter and trim                            -
##-------------------------------------------------------------------------

## the following was adapted from: https://benjjneb.github.io/dada2/tutorial.html - DADA2 Pipeline Tutorial (1.16), see the website for more information

## change the absolute path to the directory containing the fastq files
path <- here("~/Raw_data/Tick16sMiSeq2019/MiSeq_from_David/2021_Aug_cutadapt_demult_spacer_e0.1/Renamed/To_download_2021Aug20") 

fnFs <- sort(list.files(path,pattern="_16s.1.fastq.gz", full.names=TRUE)) #change the pattern to match your sample names

fnRs <- sort(list.files(path, pattern="_16s.2.fastq.gz", full.names = TRUE))

##  Extract sample names, assuming filenames have format: SAMPLENAME_###.fastq
## put "_1" to cut everything after and including "6s....)
sample.names <- sapply(strsplit(basename(fnFs), "_1"), `[`, 1)

## set up directory called filtered_trunc_trim_new for filtered sequence data
filtFs_trunc_trim <- file.path(path, "filtered_trim_new_aug29_noprim", paste0(sample.names, "_F_filt_trim_aug29_noprim.fastq.gz"))
filtRs_trunc_trim <- file.path(path, "filtered_trim_new_aug29_noprim", paste0(sample.names, "_R_filt_trim_aug29_noprim.fastq.gz"))
names(filtFs_trunc_trim) <- sample.names
names(filtRs_trunc_trim) <- sample.names

out2 <- filterAndTrim(fnFs, filtFs_trunc_trim, fnRs, filtRs_trunc_trim, maxN = 0, maxEE = c(2,2), truncQ = 2, rm.phix = TRUE, trimLeft=c(20,19), truncLen=c(240,239), compress=TRUE, multithread=TRUE)
out2

##-------------------------------------------------------------------------
##                    		Step 3 - Denoise                              -
##-------------------------------------------------------------------------

## learn error rates - 

## from time to time the self-consistency loop terminated before convergence within 10 steps, so suggest increase to 20

errF <- learnErrors(filtFs_trunc_trim, randomize = T, MAX_CONSIST = 20, verbose = T, multithread = TRUE)

errR <- learnErrors(filtRs_trunc_trim, randomize = T, MAX_CONSIST = 20, verbose = T, multithread=TRUE)

dadaFs <- dada(filtFs_trunc_trim, err=errF, pool = TRUE, verbose = TRUE, multithread=TRUE)
dadaRs <- dada(filtRs_trunc_trim, err=errR, pool = TRUE, verbose = TRUE, multithread=TRUE)


##-------------------------------------------------------------------------
##        Step 4 - Merge, size select, collapse and chimera removal       -
##-------------------------------------------------------------------------

## merge paired-end data

mergers <- mergePairs(dadaFs, filtFs_trunc_trim, dadaRs, filtRs_trunc_trim, verbose = TRUE)

## construct sequence table - this is amplicon sequence variant table
seqtab <- makeSequenceTable(mergers)

## filter the merged data for sequences of lengths ranging from 250 - 265 (V4 region is 254 bp in most, only shifting by a few basepaires in lengths)
table(nchar(getSequences(seqtab)))

## the size selection should be varied depending on the 16S primers used and what size amplicons would be considered spurious
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 250:265]
table(nchar(getSequences(seqtab2)))

seqtab2.col <-collapseNoMismatch(seqtab2, minOverlap = 20, orderBy = "abundance",   identicalOnly = FALSE, vec = TRUE, band = -1, verbose = TRUE)

## Bimera removal
seqtab.nochim <- removeBimeraDenovo(seqtab2.col, method="consensus", multithread=TRUE, verbose=TRUE)

table(nchar(getSequences(seqtab.nochim)))
table(nchar(getSequences(seqtab2.col)))

sum(seqtab.nochim)/sum(seqtab2.col)

## set pathway for RDS object output count table
saveRDS(seqtab.nochim, here("test_data/seqtab.nochim_sep5_job10.RDS")) 

##-------------------------------------------------------------------------
##                  Step 5 - Assign taxonomy and tracker                  -
##-------------------------------------------------------------------------

## accessed August 26, 2021 from - https://zenodo.org/record/4310151 (maintained on https://benjjneb.github.io/dada2/training.html)

## assign taxonomy from the rdp_train_set_18
taxa <- assignTaxonomy(seqtab.nochim, here("reference/rdp_train_set_18.fa"), tryRC=TRUE, multithread=TRUE)

## add species-level designations for 100 % matches
taxa <- addSpecies(taxa, here("reference/rdp_species_assignment_18.fa"), tryRC=TRUE)

## save taxonomy RDS object
saveRDS(taxa, here("test_data/taxa_test_sep5_job10.RDS"))

## this is a tracking tool showing read reduction throughout the four step fitlering, denoising, collapsing and chimera removal process
getN <- function(x) sum(getUniques(x))
track <- cbind(out2, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN),rowSums(seqtab2.col), rowSums(seqtab.nochim))

## If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged","collapsed", "nonchim")
rownames(track) <- sample.names
track

q()
