## ---------------------------
##
## Script name: BlastPhylo_v20220418
##
## Purpose: This is an R script used for local BLAST search against NCBI's 16S rRNA sequences (Bacteria and Archaea) and optimized phylogeny builder - see Pre-print: Paulson A.R., Huang, D., and Colautti, R.I. 2021. Evidence for symbionts and pathogens of blacklegged ticks (<i>Ixodes scapularis</i>) in an emerging Lyme disease hotspot. ##Preprint - webpage link to be added##
##
## Author: Functions by David Huang
##
## Date Created: 2021-12-06
##
##
## ---------------------------
##
## Notes: 
## 
## See the READ.ME for more information on how to get started - https://github.com/damselflywingz/tick_microbiome
## The R Project should be opened with RStudio first, and then next proceed to open and run the script.    
##
## ---------------------------

## 16S Phylogenic Analysis
# Packages :
library(Biostrings) 
library(rBLAST) 
library(genbankr) 
library(dplyr) 
#library(plyr) #this package and here are incompatible
library(ape) 
library(phangorn) 
library(seqRFLP)
library(ips)
library(rentrez)
library(here)

## Working directory :
here() # noting dplyr needs to be detached for here() to work

## NCBI BLAST software must be installed first - https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
## downloaded https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.12.0+-win64.exe Dec 1, 2021 version 2.12.0+

## set the path to the BLAST software

Sys.setenv(
  PATH = paste(
    Sys.getenv("PATH"), 
    "C:\\Program Files\\NCBI\\blast-2.12.0+\\bin", 
    sep = ";"
  )
)

## 16S Blast database https://ftp.ncbi.nlm.nih.gov/blast/db/ : 

genbank_db2 <- blast(db=(here("reference/16S_ribosomal_RNA_14Aug21V5/16S_ribosomal_RNA")))

# Samples :
samples <- readDNAStringSet(here("20220502_core_job11.fasta"))



##-------------------------------------------------------------------------
##                                Borrelia                                -
##-------------------------------------------------------------------------

# Target Sequence :
borrelia_target_ASV3 <-samples[3]
borrelia_target_ASV2 <- samples[2]

## load the functions found at the end of this script first. If the following pulls errors, try to close the RProject and re-open and run again
library(plyr)

# Blast :
blast_results <- blast_seq(borrelia_target_ASV3, genbank_db2, 32, "Borrelia Target ASV3")

detach("package:plyr", unload = TRUE)
bb1 = readDNAStringSet(here("MH781147.1.fasta"))

## add B. miyamotoi, which is a second Borrelia target found in the core microbiome that can be placed within this phylogeny
blast_results2 = add_target(blast_results,borrelia_target_ASV2,"Borrelia Target ASV2")

## generate shorter and more informative tip labels by combining Taxa and SubjectID
vec = paste(blast_results2$Taxa,blast_results$SubjectID)
vec2 = gsub("16S ribosomal RNA, ","",vec)

blast_results2$label = vec2
names_new = blast_results2$label

## the reference sequence for Borreliella burgdorferi contains 'NNN", however there are better sequences without these ambiguous bases on genbank:
## Borrellia burgdorferi isolate 15-0797 16S ribosomal RNA gene, partial sequences MH781147.1
## extracted nucleotides 399-653, based on BLASTn results, and reverse complemented
blast_results3 = add_target(blast_results2, bb1, "Borrellia burgdorferi isolate 15-0797 partial sequence. MH781147.1")

## fix the labeling with the Accession numbers
SubjectID = blast_results$SubjectID
SubjectID<-c(SubjectID,NA,NA)
vec1 = paste(blast_results3$Taxa,SubjectID)
blast_results3$label = vec1
names_new = blast_results3$label
seq_try1 = blast_results3$Seqs
names(seq_try1) = names_new
seq_try2 <- Biostrings::DNAStringSet(seq_try1)

## align the sequences
library(DECIPHER)

alignment <- AlignSeqs(DNAStringSet(seq_try2), anchor=NA)
alignment_staggered <- StaggerAlignment(alignment)
BrowseSeqs(alignment_staggered, htmlFile="20220418_borrel_stag.html", openURL = FALSE, highlight=1)

## make the NJ tree with phangorn
library(phangorn)
phang.align <- phyDat(as(alignment_staggered, "matrix"), type="DNA") #NNN is not considered

# Distance based tree :
db_tree <- optimize_db_parsimony(phang.align)

# Maximum likelihood tree :
ml_tree <- optimize_likelihood(phang.align, db_tree) # Substitution model: GTR optimize Gamma and Inv

# Bootstrap :
bs_ml_tree <- bootstrap_tree(ml_tree, 1000, 0.01, "20220418_MLBorrelia.pdf") # bootstrap label cut off is the "p" in the bootstrap_tree function, set to 65 or greater to print the node label



##-------------------------------------------------------------------------
##                                Anaplasma                               -
##-------------------------------------------------------------------------

# Target Sequence :
anaplasma_target_ASV6 <-samples[6]
AP_strains <- readDNAStringSet(here("A_phagocytophilum_strains.fasta"))
HG916766.1_Ap_ha_CSF_23 <- AP_strains[1]
HG916767.1_Ap_var_1_CSF_21 <- AP_strains[2]

library(plyr)

# Blast :
blast_results4 <- blast_seq(anaplasma_target_ASV6, genbank_db2, 11, "Anaplasma phagocytophilum Target ASV6")
detach("package:plyr", unload = TRUE)

blast_results5 = add_target(blast_results4,HG916766.1_Ap_ha_CSF_23,"HG916766.1_Ap_ha")
blast_results6 = add_target(blast_results5,HG916767.1_Ap_var_1_CSF_21,"HG916767.1_Ap_var_1")

## make the labels more informative
vec3 = paste(blast_results6$Taxa,blast_results6$SubjectID)
vec4 = gsub("16S ribosomal RNA, ","",vec3)
blast_results6$label = vec4
names_new3 = blast_results6$label
seq_try_a = blast_results6$Seqs
names(seq_try_a) = names_new3
seq_try_b <- Biostrings::DNAStringSet(seq_try_a)

# Align :
alignment_b <- AlignSeqs(DNAStringSet(seq_try_b), anchor=NA)
alignment_d_staggered <- StaggerAlignment(alignment_b)
BrowseSeqs(alignment_d_staggered, htmlFile="20220418_Annaplasma_stag.html", openURL = FALSE, highlight=1)

# Optimized phylogeny :
phang.align_d <- phyDat(as(alignment_d_staggered, "matrix"), type="DNA")
db_tree_d <- optimize_db_parsimony(phang.align_d)
ml_tree_a <- optimize_likelihood(phang.align_d, db_tree_d) # Substitution model: GTR optimize Gamma 
bs_ml_tree_a <- bootstrap_tree(ml_tree_a, 1000, 0.01, "20220418_MLAnaplasma.pdf")



##-------------------------------------------------------------------------
##                                Functions                               -
##-------------------------------------------------------------------------

## Load these first

#BLAST refseq, create DF containing BLAST results
## ap changed Taxa label truncations because part of the names getting cut off
## ap removed the "(" and ")" surrounding the Taxa label

blast_seq<-function(refseq ,blastDB,hits,targetName){
  BR<-predict(blastDB, refseq) #Blast sequence against 16S database, returns up to 500 hits
  BRhitLen<-length(BR$SubjectID) #return number of blast hits
  if(BRhitLen>hits) {
    BRhitLen=hits #Limit hit length
  }
  trun_BR<<-as.vector(BR[1:BRhitLen,])
  BRSeqs <<- DNAStringSet()
  for (i in 1:nrow(trun_BR)) {
    print(i)
    temp_gb <- readGenBank(GBAccession(trun_BR$SubjectID[i]))
    temp_dss <- temp_gb@sequence
    temp_dss@ranges@NAMES <- temp_gb@definition
    BRSeqs <<- c(BRSeqs, temp_dss)
  }
  for (i in 1:length(BRSeqs)) {
    seq_length_diff <- abs(nchar(paste(refseq)) - trun_BR[i,]$Alignment.Length)
    if (seq_length_diff > 1000) {
      trun_BR[i,]$S.end <- trun_BR[i,]$S.end + seq_length_diff
    }
  }
  tree_labels <- gsub('16S .*', '',substr(BRSeqs@ranges@NAMES,100,nchar(BRSeqs@ranges@NAMES)))
  trun_BR$Taxa<-paste0(tree_labels, substr(BRSeqs@ranges@NAMES,1,100))
  trun_BR$Seqs<-substr(paste(BRSeqs),trun_BR$S.end,trun_BR$S.start)
  targetDF<-data.frame(Taxa=targetName, Seqs=reverseComplement(refseq))
  GBDF<-rbind.fill(trun_BR,targetDF)
  GBDF<- GBDF %>% distinct(Taxa, .keep_all = TRUE)
  return(GBDF)
}

#Add more targets
add_target<-function(SeqDF,newTar,targetName){
  tar2 <-data.frame(Taxa=targetName, Seqs=reverseComplement(newTar))
  tar2DF <-merge(SeqDF,tar2, all.y = T)
  finalDF<-rbind(SeqDF,tar2DF)
}

#Align sequences 
align_seqs<-function(SeqDF){
  GbkDNA<-sapply(paste(SeqDF$Seqs),strsplit,split="")
  names(GbkDNA)<-paste(SeqDF$Taxa)
  GbkDNA<-as.DNAbin(GbkDNA)
  GbkAlign<-muscle(GbkDNA,quiet=F)
  GbkAlignTrimmed<-trimEnds(GbkAlign)
  return(GbkAlignTrimmed)
} 

#Convert alignment object to phydata type
convert_phyobject<-function(aligned_seqs){
  phy_object<-phyDat(aligned_seqs, type = "DNA")
  return(phy_object)
}

#optimize distance based methods NJ and UPGMA
optimize_db_parsimony<-function(phy_object) {
  dna_dist <- dist.ml(phy_object)
  treeNJ <- NJ(dna_dist)
  treeUPGMA <- upgma(dna_dist)
  if (parsimony(treeNJ,phy_object) < parsimony(treeUPGMA,phy_object)) {
    fit = pml(treeNJ,phy_object)
  } else {
    fit = pml(treeUPGMA,phy_object)
  }
  return(fit)
  
}

#optimize likelihood
optimize_likelihood<-function(phy_object, pml_object) {
  model_test <- modelTest(pml_object)
  model_test <- model_test %>% arrange(logLik)
  params <- strsplit(paste(model_test[nrow(model_test),]$Model),"\\+")
  model <- params[[1]][1]
  Gamma <- FALSE
  Inv <- FALSE
  if (length(params[[1]]) == 2) {
    if (params[[1]][2] == "I") {
      Inv = TRUE
    } else {
      Gamma = TRUE
    }
  } 
  if (length(params[[1]]) == 3) {
    Gamma <- TRUE
    Inv <- TRUE
  }
  print("----------------------------------------------------------------------")
  fit <- optim.pml(pml_object, model = model, optInv = Inv, optGamma = Gamma, rearrangement = "NNI")
  print("----------------------------------------------------------------------")
  print(model_test[nrow(model_test),])
  print(paste0("Gamma = ",Gamma))
  print(paste0("Inv = ", Inv))
  print("----------------------------------------------------------------------")
  print(paste0("Unoptimized loglikelihood: ",pml_object$logLik))
  print("----------------------------------------------------------------------")
  print(paste0("Optimized loglikelihood: ",fit$logLik))
  print("----------------------------------------------------------------------")
  return(fit)
}

#Build phylogenic trees
bootstrap_tree<-function(fitted_model,bs_iterations,scale_bar,out){ 
  bs <- bootstrap.pml(fitted_model, bs=bs_iterations, optNni=TRUE, multicore=F, control = pml.control(trace=0))
  pdf(out,width=13,height=10)
  plotBS(midpoint(fitted_model$tree),bs, p = 65, type="p")
  add.scale.bar(length = scale_bar, cex = 0.9, font = 2)
  dev.off()
}

## End