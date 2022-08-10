## ---------------------------
##
## Script name: 16S_decontam_analysis_Paulson_v20220401
##
## Purpose: This is an R script used to analyze 16S rRNA sequencing from the microbiome of the blacklegged tick, <i>Ixodes scapularis</i> - see Pre-print: Paulson A.R., Lougheed, S.C., Huang, D., and Colautti, R.I. 2022. Evidence for symbionts and pathogens of blacklegged ticks (<i>Ixodes scapularis</i>) in an emerging Lyme disease hotspot. ##Preprint - webpage link to be added##
##
## Author: Amber R. Paulson, Ph.D.
##
## Date Created: 2021-10-08
## Date updated: 2022-08-09
##
## Copyright (c) Amber Paulson, 2022
## Email: Amber.Rose.Paulson@gmail.com
## Github: https://github.com/damselflywingz
##
## ---------------------------
##
## Notes: 
## 
## See the READ.ME for more information on how to get started - https://github.com/damselflywingz/tick_microbiome
## !IMPORTANT - The R Project should be opened with RStudio first, and then next proceed to open and run the script.    
##
## ---------------------------

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.13")


###########################################################################
###########################################################################
###                                                                     ###
###                     SECTION 1 Decontamination                       ###
###                                                                     ###
###########################################################################
###########################################################################

##-------------------------------------------------------------------------
##                    Input RDS objects from DADA2                        -
##-------------------------------------------------------------------------

## the the ASV and taxonomy tables are provided with the RProj as RDS objects in the "data/" directory

library(phyloseq)
library(Biostrings)
library(ggplot2)
library(here) # beware plyr also has 'here' function that may need to be detached

## this is used to make relative paths based in the R project - https://www.tidyverse.org/blog/2017/12/workflow-vs-script/, https://cran.r-project.org/web/packages/here/vignettes/here.html
## this is to avoid making a bunch of absolute paths and will keep the tutorial more generic for folks that can just download the Rproject from gitHub

## here() to see where you are

## ASV table/ OTU table



seqtab.nochim_job10 <- readRDS(here("data/seqtab.nochim_sep5_job10.RDS")) #All replicates (i.e., A, B, C) PCR replicate libraries for preliminary analys
seqtab.nochim_job11 <- readRDS(here("data/seqtab.nochim_sep5_job11_AandC_only.RDS")) #A & C replicate libraries only used for analysis (i.e., B replicate omited)

## Taxonomy table

taxa_job10 <- readRDS(here("data/taxa_test_sep5_job10.RDS"))
taxa_job11 <- readRDS(here("data/taxa_test_sep5_job11_AandC_only.RDS"))

## Meta-data
samdf <- read.table(here("data/samples20.txt")) 

## need to make header = T before tidyverse can find column "Sample"
colnames(samdf) <- as.character(unlist(samdf[1,]))
samdf = samdf[-1, ]

library(tidyverse)

## this will coerce row.names=1 
samdf_try = samdf %>% 
  remove_rownames %>% column_to_rownames(var="Sample")

## make ps(sample_data$Stock_Conc) numeric
samdf_try$Stock_Conc <- as.numeric(as.character(samdf_try$Stock_Conc))

## import the phyloseq object
ps_w_b = phyloseq(otu_table(seqtab.nochim_job10, taxa_are_rows =FALSE),
               sample_data(samdf_try),
               tax_table(taxa_job10))

## move sequence from the taxa names to the refseq() slot of the phyloseq object
dna <- Biostrings::DNAStringSet(taxa_names(ps_w_b))
names(dna) <- taxa_names(ps_w_b)
ps_w_b <- merge_phyloseq(ps_w_b, dna)
taxa_names(ps_w_b) <- paste0("ASV", seq(ntaxa(ps_w_b)))


##-------------------------------------------------------------------------
##  preliminary plots to check the a, b, and c pcr replicate libraries    -
##-------------------------------------------------------------------------

ps.prelim = transform_sample_counts(ps_w_b, function(otu) otu/sum(otu))

minTotRelAbun = 0.002
x = taxa_sums(ps.prelim)
keepTaxa = taxa_names(ps.prelim)[which((x / sum(x)) > minTotRelAbun)]
pre.prunedSet = prune_taxa(keepTaxa, ps.prelim)

ps.prop.prelim = transform_sample_counts(pre.prunedSet, function(otu) otu/sum(otu))

## generate preliminary barplot for a, b, c replicates

plot <- plot_bar(ps.prop.prelim, x="Sample_full",fill = "Genus") + 
  geom_bar(aes(fill=Genus), stat="identity",width=.9) +
  facet_grid(~Location1, scales="free_x",space="free_x")

plot + theme_bw()+
  theme(axis.text.x = element_text(siz=5, angle=270)) 

ggsave(here("figs","20220401_Supplemental_relAbudn_A_B_C.pdf"), width = 9, height = 6.5)
  
library(svglite)
ggsave(here("figs","20210401_Supplemental_reAbudn_A_B_C.svg"), width = 9, height = 6.5)

## write out the ASV tax tables for Rickettsia
## identified that pcr replicate b as atypical due to batch effect compared to a, and c.

ps.prelim_rick <- subset_taxa(ps_w_b,(Genus=="Rickettsia"))

## Extract abundance matrix from the phyloseq object
prelim.asv = as(otu_table(ps.prelim_rick), "matrix")

## transpose if necessary
if(taxa_are_rows(ps.prelim_rick)){prelim.asv <- t(prelim.asv)}

## Coerce to data.frame
prelim.asv.df = as.data.frame(prelim.asv)  

write.csv(prelim.asv.df, here("figs","20220401rickettsia_A_B_C_relicatesjob10.csv"))

## ASV1, ASV88, ASV440, ASV653, ASV790, ASV1645, ASV1719 - there are 7 With A and C
## W/ 37 Rickettsia with A, B anc C: ASV1 ASV17 ASV30 ASV75 ASV299 ASV302 ASV317 ASV329 ASV370 ASV505 ASV509 ASV526 ASV547 ASV614 ASV634 ASV643 ASV680 ASV709 ASV715 ASV736 ASV741 ASV780 ASV881 ASV892 ASV1111 ASV1135 ASV1424 ASV1498 ASV1644 ASV1691 ASV1744 ASV2132 ASV2177 ASV2180 ASV2243 ASV2395 ASV2459

library(reltools)
save_fasta(ps.prelim_rick, file = "20220401_A_B_C_replicates_rickettsia_job10.fasta", rank = "Genus")  

## based on preliminary analysis, the 'b' replicate samples were identified to contain increased Rickettsia diversity that were not observed in replicate a and c samples
## the increased diversity of Rickettsia was also evident in the PCR control sample for 'b', indicating a batch pcr effect associated with replicate b

## b was removed and the dada2 repeated on the 'a' and 'c' replicates (i.e, 'job11')

## import the phyloseq object
ps <- phyloseq(otu_table(seqtab.nochim_job11, taxa_are_rows =FALSE),
               sample_data(samdf_try),
               tax_table(taxa_job11))
ps

## job11 = 2332

## move sequence from the taxa names to the refseq() slot of the phyloseq object
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))


##-------------------------------------------------------------------------
##                            Decontaminate                               -
##-------------------------------------------------------------------------

## decontamination is completed using the decontam package following - https://benjjneb.github.io/decontam/vignettes/decontam_intro.html
library(decontam)

## tell decontam package which samples are negative control
sample_data(ps)$is.neg<-sample_data(ps)$Whole_or_dissected == "Not_applicable"
## sample_data(ps_w_b)$is.neg<-sample_data(ps_w_b)$Whole_or_dissected == "Not_applicable"

## Run the next steps of the analysis on the "ps" data set tht did not include the B replicate at all.

## separate the replicates
rep_A = subset_samples(ps,Replicate=="a")
rep_C = subset_samples(ps,Replicate %in% c("c","d"))

## repA
repA_contamdf.prev <- isContaminant(rep_A, method="prevalence", neg=sample_data(rep_A)$is.neg)
table(repA_contamdf.prev$contaminant) #7

tax_2repA <- as(tax_table(rep_A), "matrix")
contam_ASV_2repA = paste0("ASV",which(repA_contamdf.prev$contaminant))
tax_contam_2repA <- subset(tax_2repA,rownames(tax_2repA) %in% contam_ASV_2repA)
tax_contam_2repA

tax_A <- as(tax_table(rep_A), "matrix")
contam_ASV_repA = paste0("ASV",which(repA_contamdf.prev$contaminant))
tax_contam_repA <- subset(tax_A,rownames(tax_A) %in% contam_ASV_repA)
tax_contam_repA

## repC
repC_contamdf.prev <- isContaminant(rep_C, method="prevalence", neg=sample_data(rep_C)$is.neg)
table(repC_contamdf.prev$contaminant) #26

tax_2repC <- as(tax_table(rep_C), "matrix")
contam_ASV_2repC = paste0("ASV",which(repC_contamdf.prev$contaminant))
tax_contam_2repC <- subset(tax_2repC,rownames(tax_2repC) %in% contam_ASV_2repC)
tax_contam_2repC

tax_C <- as(tax_table(rep_C), "matrix")
contam_ASV_repC = paste0("ASV",which(repC_contamdf.prev$contaminant))
tax_contam_repC <- subset(tax_C,rownames(tax_C) %in% contam_ASV_repC)
tax_contam_repC

rep_A_decon <- prune_taxa(!repA_contamdf.prev$contaminant,rep_A)

rep_C_decon <- prune_taxa(!repC_contamdf.prev$contaminant,rep_C)

## re-combine the decontaminated replicate batches

ps_2 = merge_phyloseq(rep_A_decon,rep_C_decon)

ps_2 #2331

## now the "frequency" method of decontam is deployed
ps2_contamdf.freq <- isContaminant(ps_2, method="frequency",conc="Stock_Conc")
table(ps2_contamdf.freq$contaminant) #272

tax_2 <- as(tax_table(ps_2), "matrix")
contam_ASV_2 = paste0("ASV",which(ps2_contamdf.freq$contaminant))
tax_contam_2 <- subset(tax_2,rownames(tax_2) %in% contam_ASV_2)
dim(tax_contam_2)

ps2_decontam <- prune_taxa(!ps2_contamdf.freq$contaminant,ps_2)
ps2_decontam #2059

total_contam <- rbind(tax_contam_repA, tax_contam_repC, tax_contam_2)
write.csv(total_contam, here("figs/contaminating_ASVs_freq_20220401.csv"))

## extract the refseq as sequence data before merge
dna_1 <- Biostrings::DNAStringSet(refseq(ps2_decontam))

## https://github.com/benjjneb/dada2/issues/745
library(dplyr)

## the samples with only one representative library are getting zeroed out, so extact those from the phyloseq object
## then run the function on those samples with 2 replicates, then rejoin it togeter

head(sample_data(ps2_decontam)) # can subset all but, MiSeqID == c("B5","D3","E2","F2","F5","H5")
ps2_decontam_dupl = subset_samples(ps2_decontam,Sample_AC!="232_sg")
ps2_decontam_dupl = subset_samples(ps2_decontam_dupl,Sample_AC!="42_sg")
ps2_decontam_dupl = subset_samples(ps2_decontam_dupl,Sample_AC!="24_g")
ps2_decontam_dupl = subset_samples(ps2_decontam_dupl,Sample_AC!="24_sg")
ps2_decontam_dupl = subset_samples(ps2_decontam_dupl,Sample_AC!="260_sg")
ps2_decontam_dupl = subset_samples(ps2_decontam_dupl,Sample_AC!="72_sg")
ps2_decontam_dupl = subset_samples(ps2_decontam_dupl,Sample_AC!="44_g")
ps2_decontam_dupl = subset_samples(ps2_decontam_dupl,Sample_AC!="NB")
ps2_decontam_dupl = subset_samples(ps2_decontam_dupl,Sample_AC!="NP")
ps2_decontam_dupl

ps2_decontam_nodupl1 = subset_samples(ps2_decontam, Sample_AC=="232_sg")
ps2_decontam_nodupl2 = subset_samples(ps2_decontam, Sample_AC=="42_sg")
ps2_decontam_nodupl3 = subset_samples(ps2_decontam, Sample_AC=="24_g")
ps2_decontam_nodupl4 = subset_samples(ps2_decontam, Sample_AC=="24_sg")
ps2_decontam_nodupl5 = subset_samples(ps2_decontam, Sample_AC=="260_sg")
ps2_decontam_nodupl6 = subset_samples(ps2_decontam, Sample_AC=="72_sg")
ps2_decontam_nodupl9 = subset_samples(ps2_decontam, Sample_AC=="44_g")
ps2_decontam_nodupl7 = subset_samples(ps2_decontam, Sample_AC=="NB")
ps2_decontam_nodupl8 = subset_samples(ps2_decontam, Sample_AC=="NP")

## following adapted from https://github.com/benjjneb/dada2/issues/745 takes any samples represented by two libraries and only keeps ASVs that occur in both and zeros everything else.

library(dplyr)

## Merge the samples that have the same Specimen by adding reads
ps.merged_ps2_decontam_dupl <- ps2_decontam_dupl %>%
  merge_samples(group = "Sample_AC")
ps.merged_ps2_decontam_dupl_ori = ps.merged_ps2_decontam_dupl

## Create a matrix of 0s and 1s indicating whether the taxon count should be allowed, or should be set to 0. 
## changed to == 2, because there are 2 replicates (b replicate not included in analysis)
ps.merged.ok_dupl <- ps2_decontam_dupl %>%
  transform_sample_counts(function (x) (x > 0) * 1) %>%
  merge_samples(group = "Sample_AC") %>%
  transform_sample_counts(function (x) (x == 2) * 1)
ps.merged.ok_dupl

## Multiply the counts by the 0-1 matrix
newotu <- otu_table(ps.merged_ps2_decontam_dupl) * otu_table(ps.merged.ok_dupl)

## Turn newotu back into an otu_table and update the otu_table of ps.merged object
otu_table(ps.merged_ps2_decontam_dupl) <- otu_table(newotu, taxa_are_rows = FALSE)

## add the single samples back in
ps.merged_ps2_decontam_3 = merge_phyloseq(ps.merged_ps2_decontam_dupl, ps2_decontam_nodupl1,ps2_decontam_nodupl2, ps2_decontam_nodupl3, ps2_decontam_nodupl4, ps2_decontam_nodupl5, ps2_decontam_nodupl6, ps2_decontam_nodupl9, ps2_decontam_nodupl7, ps2_decontam_nodupl8)

## re-import the refseq
ps.merged_ps2_decontam_4 = merge_phyloseq(ps.merged_ps2_decontam_3,dna_1)

## check the results
head(otu_table(ps.merged_ps2_decontam_4)[, 1:5])
head(otu_table(ps.merged_ps2_decontam_dupl_ori)[,1:5])

## the single replicate samples are returned
(sample_sums(ps.merged_ps2_decontam_4)) 





###########################################################################
###########################################################################
###                                                                     ###
###                 SECTION 2 Analysis with phyloseq                    ###
###                                                                     ###
###########################################################################
###########################################################################


##-------------------------------------------------------------------------
##                        Preliminary ordination                          -
##-------------------------------------------------------------------------

## return the metadata following the merge by adding back the necessary variables to the sample_data
new_merge_data <- read.table(here("data","2022.02.15_mergeAC_datasheet.txt"))

## add the variables to the phyloseq object sample_data
AC = new_merge_data$Sample_AC
Location = new_merge_data$Location1
Tissue = new_merge_data$Type
location = new_merge_data$Location
Type = new_merge_data$Whole_or_dissected

sample_data(ps.merged_ps2_decontam_4)$Sample_AC = AC
sample_data(ps.merged_ps2_decontam_4)$Location1 = Location
sample_data(ps.merged_ps2_decontam_4)$Tissue = Tissue
sample_data(ps.merged_ps2_decontam_4)$Location=location
sample_data(ps.merged_ps2_decontam_4)$Type=Type

## remove samples with less than 5000 sequences so we can proceed
sample_sums(ps.merged_ps2_decontam_4)

## export tax table
tax.clean = data.frame(tax_table(ps.merged_ps2_decontam_4))
tax.clean[tax.clean == "-"] = NA ## custom for rdp as there are some problematic entries "-"

tax.clean[is.na(tax.clean)] <- "" ## then replace all NA with blank

## the next part is adapted from https://github.com/joey711/phyloseq/issues/850, which tidies up the taxonomy table.

## change columns to characters
for (i in 1:7){ tax.clean[,i] <- as.character(tax.clean[,i])}

tax.clean[is.na(tax.clean)] <- ""

for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,2] == ""){
    kingdom <- paste(tax.clean[i,1]," (kingdom)",  sep = "")
    tax.clean[i, 2:7] <- kingdom
  } else if (tax.clean[i,3] == ""){
    phylum <- paste(tax.clean[i,2]," (phylum)",  sep = "")
    tax.clean[i, 3:7] <- phylum
  } else if (tax.clean[i,4] == ""){
    class <- paste(tax.clean[i,3]," (class)",  sep = "")
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i,5] == ""){
    order <- paste(tax.clean[i,4], " (order)",  sep = "")
    tax.clean[i, 5:7] <- order
  } else if (tax.clean[i,6] == ""){
    family <- paste(tax.clean[i,5]," (family)",  sep = "")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] != ""){
    tax.clean$Species[i] <- paste(tax.clean$Genus[i]," ",tax.clean$Species[i], sep = "")
  } else if (tax.clean[i,7] == ""){
    tax.clean$Species[i] <- paste(tax.clean$Genus[i], " sp.", sep = "")
  }
}

## return to phyloseq object
tax_table(ps.merged_ps2_decontam_4) <- as.matrix(tax.clean)
head(tax_table(ps.merged_ps2_decontam_4), n= 10)

ps_ <- ps.merged_ps2_decontam_4

## typical reproducible microbiome workflows will removal suspected non-bacterial taxa not assigned at phylum-level
ps0 <- subset_taxa(ps_, !is.na(Phylum) & !Phylum %in% c("", "Bacteria (kingdom)"))
                   
head(tax_table(ps0), n= 10) ## ASV7 is gone, Babesia apicoplast interference - this was not consistently recovered across replicates so zero sample_sums

## ps0 job11 18942

ps.Bacteria <- subset_taxa(ps_, (Phylum == "Bacteria (kingdom)"))
head(otu_table(ps.Bacteria)[, 1:5])
sample_sums(ps.Bacteria)
ps.Bacteria #155

#library(reltools)
ps.ASV7 <- prune_taxa("ASV7", ps_)
sample_sums(ps.ASV7)
save_fasta(ps.ASV7, file = "20220407_ASV7.fasta", rank = "Phylum") # F5_a and H5_a

## blasted the sequence against nt database, top hits are for Babeisa sp. identified from sheep in China (Unpublished; MH992225.1, MH992224.1,  KX881914.1; all 94 % - 248/264), Babesia orientalis, Babeisa bigemina, Babesia motasi, and Babesia sp.

m = taxa_sums(ps0)
taxa_names(ps0)[which(m <= 4)] ## 810 ASVs 

ps0. = prune_taxa(taxa_names(ps0)[which(m > 4)], ps0) ## only keep taxa with taxa_sum greater than 4
ps0. # 1094

## remove non-bacterial ASVs

ps0.1 <- subset_taxa(ps0., (Kingdom!="Archaea") | is.na(Kingdom))
ps0.1
## 1088

ps0.2 <- subset_taxa(ps0.1, (Phylum!="Cyanobacteria/Chloroplast") | is.na(Phylum))
ps0.2
## 1062

ps0.3 <- subset_taxa(ps0.2, (Kingdom!="") | is.na(Kingdom))
ps0.3

## 1048

ps0.4 <- subset_taxa(ps0.3, (Kingdom!="Eukaryota") | is.na(Kingdom))
ps0.4 #1042

## ps0.5 <- subset_taxa(ps0.4, (Phylum!="Mitochondria") | is.na(Phylum)) # no Mitichondria present

Frans <- subset_taxa(ps0.2, (Genus == "Francisella")) ## ASV216 (endosymbiont), ASV297 (inquiline) 
tax_table(Frans)
sample_sums(Frans) ## 72_g, H2_c (trace), F5_a, H5_a
save_fasta(Frans, file = "20220406_Francisella.fasta", rank = "Species")  

## alignment by BLAST 98 % (254/258) identities
## see Table S2 and Table S5* in Rynkiewicz et al. 2015 (https://onlinelibrary.wiley.com/doi/10.1111/mec.13187)
## some low prevelance of FLE detected here in IS ticks 
## lower relative abundance Francisella is an inquiline.

Rick_A_C <- subset_taxa(ps0.2, (Genus == "Rickettsia")) # 6 taxa
sample_sums(Rick_A_C)

save_fasta(Rick_A_C, file = "20220407_Rickettsia_6ASVs.fasta", rank = "Species") 




##-------------------------------------------------------------------------
##                              Rarefaction                               -
##-------------------------------------------------------------------------

## remove all negatives
ps_noneg <- prune_samples(sample_sums(ps0.4) >= 7500, ps0.4)
sample_sums(ps_noneg)

## rarefaction wrapper ggrare of the ranacapa package, based on vegan rarecurve 
library(vegan)

#if (!requireNamespace("devtools", quietly = TRUE))
 # install.packages('devtools')

library(devtools)

#devtools::install_github("gauravsk/ranacapa")
library(ranacapa)

rare_graph <- ggrare(ps_noneg, step = 50, label = "Sample_AC", se = FALSE)
rare_graph + theme_bw() + guides(fill="none", color="none")

## to calculate the number of ASVs detected in each sample followed instructions from https://github.com/joey711/phyloseq/issues/337

## Create a factor corresponding to the Genera
genfac = factor(tax_table(ps_noneg)[, "Species"])

## Tabulate the counts for each genera in each sample
gentab = apply(otu_table(ps_noneg), MARGIN = 1, function(x) {
  tapply(x, INDEX = genfac, FUN = sum, na.rm = TRUE, simplify = TRUE)
})

## your threshold, in your case, 1.
observationThreshold = 1
apply(gentab > observationThreshold, 2, sum)

ggsave(here("figs","20220407Rarefaction_supplement.pdf"),  width = 9, height = 6.5)
ggsave(here("figs","202204073Rarefaction_supplement.svg"),  width = 9, height = 6.5)


##-------------------------------------------------------------------------
##                  Characterization of core microbiome                   -
##-------------------------------------------------------------------------

ps0.4_prop = transform_sample_counts(ps0.4, function(otu) otu/sum(otu)) ## apply total sum scaling normalization (i.e., relative abundance)

## characterization of chore microbiome (greater than 0.1 % relative abundance)
## https://github.com/joey711/phyloseq/issues/694

minTotRelAbun = 0.001
x = taxa_sums(ps0.4_prop)
keepTaxa = taxa_names(ps0.4_prop)[which((x / sum(x)) > minTotRelAbun)]
prunedSet = prune_taxa(keepTaxa, ps0.4_prop) ## 36 remain now  

prune.pos = transform_sample_counts(prunedSet, function(otu) otu/sum(otu))

## define the number of colors based on number of ASVs identified in the ps object
genusList = unique(tax_table(prune.pos)[,"Genus"])
speciesList = unique(tax_table(prune.pos)[,"Species"])

x = dim(genusList)
nb.cols.genus=x[1]

xy = dim(speciesList)
nb.cols.sp=xy[1]

library(randomcoloR)

palette <- distinctColorPalette(nb.cols.sp)

## Generate supplmentary figure showing total abundance of negative control and sample libraries
ps0.4_pruned = prune_taxa(keepTaxa, ps0.4)
tax_table(ps0.4_pruned) <- as.matrix(tax.clean)
sample_sums(ps0.4_pruned)

library(scales)

## the thin black outline is obscuring the top of the barplot, work around adapted from: https://stackoverflow.com/questions/52747802/how-to-remove-very-thin-bar-plot-outline-border
plot_bar_2 <-  function (physeq, x = "Sample_AC", y = "Abundance", fill = NULL, title = NULL, facet_grid = NULL, border_color = NA) 
{
  mdf = psmelt(physeq)
  p = ggplot(mdf, aes_string(x = x, y = y, fill = fill))
  p = p + geom_bar(stat = "identity", position = "stack",  color = border_color)
  p = p + theme(axis.text.x = element_text(angle = -90, hjust = 0))
  if (!is.null(facet_grid)) {
    p <- p + facet_grid(facet_grid)
  }
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  return(p)
}

p <- plot_bar_2(ps0.4_pruned, fill = "Species") +
  facet_grid(~Location1, scales="free_x",space="free_x") +
  scale_y_continuous(labels =comma)

p + geom_bar(stat = "identity") + scale_fill_manual(values = palette)+
theme_bw()+
  xlab("Library")+
  theme(axis.text.x = element_text(siz=8, angle=270))

ggsave(here("figs","20220408_Supplemental_totAbudn_bar.pdf"), width = 9, height = 6.5)
ggsave(here("figs","20220408_Supplemental_totAbudn_bar.svg"), width = 9, height = 6.5)

## repeat on the noneg object
ps_noneg_prop = transform_sample_counts(ps_noneg, function(otu) otu/sum(otu))

minTotRelAbun = 0.001
x = taxa_sums(ps_noneg_prop)
keepTaxa = taxa_names(ps_noneg_prop)[which((x / sum(x)) > minTotRelAbun)]
prunedSet_noneg = prune_taxa(keepTaxa, ps_noneg_prop) ## 25 remain now  

prune.pos_noneg = transform_sample_counts(prunedSet_noneg, function(otu) otu/sum(otu))

## define the number of colors based on number of ASVs identified in the ps object
genusList = unique(tax_table(prune.pos_noneg)[,"Genus"])
speciesList = unique(tax_table(prune.pos_noneg)[,"Species"])

x = dim(genusList)
nb.cols.genus=x[1]

xy = dim(speciesList)
nb.cols.sp=xy[1]

ps.filt2_pruned_noneg = prune_taxa(keepTaxa, ps_noneg)

library(reltools)
save_fasta(ps.filt2_pruned_noneg, file = "20220502_core_job11.fasta", rank = "Species")

## subset without control samples - see total abundance plot below, the negative PCR control had less sequences compared to sample libraries and the negative bead control had unique ASVs not identified in other libraries

mycolors_outline <- c("#6d6d6e")
palette <- distinctColorPalette(nb.cols.sp)
## plot the relative abundance of core microbiome without controls

p <- plot_bar_2(prune.pos_noneg, fill = "Species") +
  facet_grid(~Location1, scales="free_x",space="free_x") +
  scale_y_continuous(labels =comma) + 
  scale_fill_manual(values = palette)+
  theme_bw()

p1 = p + geom_bar(aes(fill=Species,color=""), size = .1,stat="identity",position = "stack") +
  scale_color_manual(values=mycolors_outline)+
  theme_bw()

## increase text size for publication
p1 + labs(y = "Relative abundance", x = "V4 16S rDNA library")+
  theme(axis.text.x = element_text(siz=11, angle=270,vjust=0),
        axis.text.y.left = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.title = element_text(size=12),
        legend.text = element_text(size = 12)) 

ggsave(here("figs","20220408_relAbudn_bar.pdf"), width = 9, height = 6.5)
ggsave(here("figs","20220408_relAbudn_bar.svg"), width = 9, height = 6.5)


##-------------------------------------------------------------------------
##            violin and boxplot plot of ASVs of concern                  -
##-------------------------------------------------------------------------

## generate the otu_table with relative abundance of core microbiome and then make 

top3 <- names(sort(taxa_sums(prune.pos_noneg), decreasing=TRUE))[1:3]
Ana <- subset_taxa(prune.pos_noneg, (Genus == "Anaplasma"))
tax_table(Ana)
ps.top3 <- prune_taxa(top3, prune.pos_noneg)
ps.top3_A <- merge_phyloseq(ps.top3, Ana)
ps.top3_A

## subset without 'rem' and control samples - the 'rem' is removed because only 3 libraries and not informative here
ps.top3_sub = subset_samples(ps.top3,Tissue!="rem")
ps.top3_A_sub = subset_samples(ps.top3_A,Tissue!="rem")
Top3_data_A = as(otu_table(ps.top3_A_sub), "matrix")

## transpose
if(taxa_are_rows(ps.top3_A_sub)){Top3_data_A <- t(Top3_data_A)}

## coerce to a dataframe3
Top3_dfA = as.data.frame(Top3_data_A)

detach("package:plyr", unload = TRUE)
write.csv(Top3_dfA, here("top3_Ana_otu_relabund_20220408.csv"))

library(tidyverse)

## this will coerce row.names=1 
Top3_dfA = cbind(rownames(Top3_dfA), data.frame(Top3_dfA, row.names=NULL))

#install.packages("~/R/R-4.1.1/library/plyr_1.8.6.zip", repos = NULL, type = "win.binary")
library(plyr)

Top3_dfA_ = rename(Top3_dfA, c("rownames(Top3_dfA)"="Sample"))
data.frame(append(Top3_dfA_, list(Pathogen=c(rep(c("ASV1"), times=33))), after=match("ASV8", names(Top3_dfA_))))
Sample_repA = rep(c(Top3_dfA_$Sample), times = 4)

## here need to replace the single sample values see: https://statisticsglobe.com/r-replace-value-of-data-frame-variable-dplyr-package
sam1 = replace(Sample_rep, Sample_rep=="B5_a","232_sg")
sam2 = replace(sam1, sam1=="E2_a","24_g")
sam3 = replace(sam2, sam2=="F2_a","24_sg")
sam4 = replace(sam3, sam3=="F5_a","260_sg")
sam5 = replace(sam4, sam4=="H5_a","72_sg")
sam6 = replace(sam5, sam5=="H2_c","44_g")
Sample_rep = sam6

sam1 = replace(Sample_repA, Sample_repA=="B5_a","232_sg")
sam2 = replace(sam1, sam1=="E2_a","24_g")
sam3 = replace(sam2, sam2=="F2_a","24_sg")
sam4 = replace(sam3, sam3=="F5_a","260_sg")
sam5 = replace(sam4, sam4=="H5_a","72_sg")
sam6 = replace(sam5, sam5=="H2_c","44_g")
Sample_repA = sam6

Pathogen_repA = rep(c("ASV1","ASV2","ASV3","ASV8"), each=33)

Rel_combA = c(Top3_dfA_$ASV1, Top3_dfA_$ASV2, Top3_dfA_$ASV3, Top3_dfA_$ASV8)

Top3_mergeA = cbind(data.frame(Sample_repA, Pathogen_repA,Rel_combA))

library(dplyr)
library(tidyr)

Top3_merge2 = Top3_merge %>%
  separate(Sample_rep, c("Sample_rep","type"), "_")

Top3_merge3 = Top3_merge2 %>% 
  mutate_at(vars(Rel_comb), as.numeric)

Top3_merge2A = Top3_mergeA %>%
  separate(Sample_repA, c("Sample_rep","type"), "_")

Top3_merge3A = Top3_merge2A %>% 
  mutate_at(vars(Rel_combA), as.numeric)

#install.packages("hrbrthemes")
library(hrbrthemes)

library(viridis)

## violin plot

ggplot(Top3_merge3A, aes(x=Pathogen_repA, y=Rel_combA, fill=Pathogen_repA)) +
  geom_violin(bw=.03, size = .001) +
  geom_point(size=.5)+
  scale_fill_viridis(discrete = TRUE, alpha=0.5, option="D")+
  facet_grid(~type, scales="free_x",space="free_x")+
  xlab("Pathogens")+
  ylab("Relative abundace")+
  theme_bw()

detach("package:plyr", unload = TRUE) #must detach this package to use here

ggsave(here("figs", "20220408_pathogen_violin_relabund.pdf"), width = 9, height = 6.5)
ggsave(here("figs", "20220408_pathogen_violin_relabund.svg"), width = 9, height = 6.5)

## boxplot

ggplot(Top3_merge3A, aes(x=Pathogen_repA, y=Rel_combA, fill=Pathogen_repA)) +
  geom_boxplot(outlier.shape = NA) + # outliers are hidden, NOT removed
  geom_point(size=1)+
  scale_fill_viridis(discrete = TRUE, alpha=0.5, option="D")+
  facet_grid(~type, scales="free_x",space="free_x")+
  xlab("Pathogens")+
  ylab("Relative abundace")+
  theme_bw()

ggsave(here("figs", "20220408_pathogen_Ana_box_relabund.pdf"), width = 9, height = 6.5)
ggsave(here("figs", "20220408_pathogen_Ana_box_relabund.svg"), width = 9, height = 6.5)


##-------------------------------------------------------------------------
##                UniFrac ordination and outlier removal                  -
##-------------------------------------------------------------------------

## working with total sum scaled (i.e., relative abundance) data
ps_noneg_prop

## 1042 and 36 samples 

## sort top 500 most abundant ASVs
top500_rel <- names(sort(taxa_sums(ps_noneg_prop), decreasing=TRUE))[1:500]
ps.top500_rel <- prune_taxa(top500_rel, ps_noneg_prop) # keep top 500 most relatively abundant ASVs, then prune from the whole count phyloseq object for Unifrac

## export a fasta file from the ps.top500_rel phyloseq object
ps.filt2_pruned500_rel = prune_taxa(top500_rel, ps_noneg)

## now make phylogeny
seqs <- refseq(ps.filt2_pruned500_rel)

## BiocManager::install("DECIPHER")
library(DECIPHER)
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)
alignment_staggered <- StaggerAlignment(alignment)

## make the NJ tree with phangorn
library(phangorn)
phang.align <- phyDat(as(alignment_staggered, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)
fitGTR <- update(fit, k=4, inv=0.2)

## this step is the computational bottle-neck 
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))

## add the phylogeny to the object
physeq1_rel = merge_phyloseq(ps.top500_rel,fitGTR$tree)

## create the weighted (considers abundance) and unweighted (presence/absence) UniFrac distances
ordu1 = ordinate(physeq1_rel, "PCoA", "unifrac", weighted=FALSE)
ordu2 = ordinate(physeq1_rel, "PCoA", "unifrac", weighted=TRUE)

plot_uw = plot_ordination(physeq1_rel, ordu1, type="samples",shape="Type",color="Location", title="Unweighted UniFrac top 500")+
  geom_point(aes(fill=Location))+
  scale_shape_manual(values = c(16, 17))+
  scale_color_manual(values = c("#00AFBB","#994bc2","#ED4F50"))+
  geom_text(mapping = aes(label=Sample_AC),size=4,nudge_x=0.01,nudge_y=-0.02)+
  theme_classic() +  #changes theme, removes grey background
  theme(legend.text = element_text(size = 12))+
  theme(axis.text.y.left = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.title = element_text(size=12))+
  guides(shape = guide_legend(override.aes = list(size=4)))+
  guides(fill = guide_legend(override.aes = list( size=4 )))

plot_uw2 = plot_uw + geom_point(aes(shape = Type, color= Location, fill = Location),size=4)
plot_uw2

# install.packages("ggforce")
library(ggforce)

plot_uw3 = plot_uw2 +
  ggforce::geom_mark_ellipse(show.legend=F)+
  coord_equal()

plot_uw3

ggsave(here("figs","20220412_AC_merged_unweighted_unifrac_top500.pdf"), width = 9, height = 6.5)
ggsave(here("figs","20220412_AC_merged_unweighted_unifrac_top500.svg"), width = 9, height = 6.5)

plot_w = plot_ordination(physeq1_rel, ordu2, type="samples",shape="Type",color="Location", title="Weighted UniFrac")+
  geom_point(aes(fill=Location))+
  scale_shape_manual(values = c(16, 17))+
  scale_color_manual(values = c("#00AFBB","#994bc2","#ED4F50"))+
  geom_text(mapping = aes(label=Sample_AC),size=4,nudge_x=0.01,nudge_y=-0.02)+
  theme_classic() +  #changes theme, removes grey background
  theme(legend.text = element_text(size = 12))+
  theme(axis.text.y.left = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.title = element_text(size=12))+
  guides(shape = guide_legend(override.aes = list(size=4)))+
  guides(fill = guide_legend(override.aes = list( size=4 )))

plot_w1 = plot_w + geom_point(aes(shape = Type, color= Location, fill = Location),size=4)

plot_w1

ggsave(here("figs","20220411_AC_merged_weighted_unifrac_top500.pdf"), width = 9, height = 6.5)
ggsave(here("figs","20220411_AC_merged_weighted_unifrac_top500.svg"), width = 9, height = 6.5)

##-------------------------------------------------------------------------
##                Agglomeration and Unifrac Ordination                    -
##-------------------------------------------------------------------------

## agglomerate at Species level - http://web.stanford.edu/class/bios221/MicrobiomeWorkflowII.html
ps.filt2_merged_glom2 = tax_glom(ps_noneg, "Species", NArm = TRUE) 
sample_data(ps.filt2_merged_glom2)
## 419 taxa and 36 samples

## Unifrac on the agglomerated ASV dataset
taxa_names <- names(sort(taxa_sums(ps.filt2_merged_glom2), decreasing=TRUE))

taxa_names_w_seq <- prune_taxa(taxa_names, ps0.4) ## keep taxa to retrieve the sequence data

## now make phylogeny
seqs_x <- refseq(taxa_names_w_seq) ## 419

## align the sequences with DECIPHER

alignment_x <- AlignSeqs(DNAStringSet(seqs_x), anchor=NA)
alignment_staggered_x <- StaggerAlignment(alignment_x)

## make the NJ tree with phangorn

phang.align <- phyDat(as(alignment_staggered_x, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)

fitGTR_all <- update(fit, k=4, inv=0.2)

## this step is the computational bottle-neck
fitGTR_all <- optim.pml(fitGTR_all, model="GTR", optInv=TRUE, optGamma=TRUE,
                        rearrangement = "stochastic", control = pml.control(trace = 0))

## Once you have created a tree object it needs to be merged back with the original phyloseq object in the tree slot (tree tip labels need to match the phyloseq object)

## merge the tree with the phyloseq object
physeq2_red = merge_phyloseq(ps.filt2_merged_glom2,fitGTR_all$tree) #hopefully can merge to the glommed/merged now ps.filt5_merged_glom
physeq2_red
## 419 taxa and 36 samples

## now that the phylogenetic tree is incorporated into the physeq1 object we can use phyloseq to do the Unifrac ordination
## scripts from: https://joey711.github.io/phyloseq/plot_ordination-examples.html

## ordinate using unifrac distances, weighted considers abundances, unweighted is presence/absence but no abundance)
ordu1_red = ordinate(physeq2_red, "PCoA", "unifrac", weighted=FALSE)
ordu2_red = ordinate(physeq2_red, "PCoA", "unifrac", weighted=TRUE)

## plot the unweighted and weighted Unifrac ordinations

qq =plot_ordination(physeq2_red, ordu1_red, type="samples",shape="Type",color="Location", title="Unweighted UniFrac agglomerated")+
  geom_point(aes(fill=Location))+
  scale_shape_manual(values = c(16, 17))+
  scale_color_manual(values = c("#00AFBB","#994bc2","#ED4F50"))+
  geom_text(mapping = aes(label=Sample_AC),size=4,nudge_x=0.01,nudge_y=-0.02)+
  theme_classic() +  #changes theme, removes grey background
  theme(legend.text = element_text(size = 12))+
  theme(axis.text.y.left = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.title = element_text(size=12))+
  guides(shape = guide_legend(override.aes = list(size=4)))+
  guides(fill = guide_legend(override.aes = list( size=4 )))

qq1 = qq+geom_point(aes(shape = Type, color= Location, fill = Location),size=4)

#install.packages("ggforce")
library(ggforce)

qq2 = qq1 +
  ggforce::geom_mark_ellipse(show.legend=F)+
  coord_equal()
qq2

ggsave(here("figs","20220411_AC_merged_unweighted_unifrac_filt2_reduced_aglommerateda.pdf"))
ggsave(here("figs","20220411_AC_merged_unweighted_unifrac_filt2_reduced_aglommerateda.svg"))

## plot the weighted UniFrac
qqq =plot_ordination(physeq2_red, ordu2_red, type="samples",color="Location",shape="Type", title="Weighted UniFrac")+
  geom_point(aes(fill=Location))+
  scale_shape_manual(values = c(16, 17))+
  scale_color_manual(values = c("#00AFBB","#994bc2","#ED4F50"))+
  geom_text(mapping = aes(label=Sample_AC),size=4,nudge_x=0.01,nudge_y=-0.02)+
  theme_classic() +                                                      #changes theme, removes grey background
  theme(legend.text = element_text(size = 12))+
  theme(axis.text.y.left = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.title = element_text(size=12))+
  guides(shape = guide_legend(override.aes = list(size=4)))+
  guides(fill = guide_legend(override.aes = list( size=4 )))

qqq+geom_point(aes(shape = Type, color= Location, fill = Location),size=4)

ggsave(here("figs","20220411_AC_merged_weighted_unifrac_filt2_reduced.pdf"), width = 9, height = 6.5)
ggsave(here("figs","20220411_AC_merged_weighted_unifrac_filt2_reduced.svg"), width = 9, height = 6.5)


##-------------------------------------------------------------------------
##          Compare distiances within and between groups                  -
##-------------------------------------------------------------------------

## inspiration from https://msystems.asm.org/content/5/3/e00292-20.abstract#sec-2 fig 1D

u_u = phyloseq::distance(physeq2_red,m="unifrac", type=c("samples")) #unifrac is unweighted, wunifrac is weighted

library(reshape2)
u_u.m = melt(as.matrix(u_u))

## following adapted from - https://github.com/joey711/phyloseq/issues/714

## remove self-comparisons
u_u.m = u_u.m %>%
  filter(as.character(Var1) != as.character(Var2)) %>%
  mutate_if(is.factor,as.character)

## order by value
u_u.m2 = u_u.m[order(u_u.m$value),]

## remove every other row (these are duplicates from vice-versa comparisons) #https://statisticsglobe.com/select-odd-even-rows-columns-from-data-frame-r

row_odd <- seq_len(nrow(u_u.m2)) %% 2
row_odd

u_u.m3 = u_u.m2[row_odd ==1, ]

## need to pull from sample data vectors "Sample_AC", and "Type"
test.1 = sample_data(physeq2_red)$Sample_AC

AC = test.1
AC1 = replace(AC, AC=="232_sg","B5_a")
AC2 = replace(AC1, AC1=="24_g","E2_a")
AC3 = replace(AC2, AC2=="24_sg","F2_a")
AC4 = replace(AC3, AC3=="260_sg","F5_a")
AC5 = replace(AC4, AC4=="72_sg","H5_a")
AC6 = replace(AC5, AC5=="44_g","H2_c")
AC = AC6

test.2 = sample_data(physeq2_red)$Type
Type = test.2

sd = melt(data.frame(AC,Type))

## combined distances with sample data
## mutating joins adds coluns y to x, matching rows based on the keys (i.e., by "Var1") - https://dplyr.tidyverse.org/reference/mutate-joins.html
colnames(sd) = c("Var1", "Type1")

u.sd1 = left_join(u_u.m3, sd, by = "Var1")
colnames(sd) = c("Var2", "Type2")

u.sd2 = left_join(u.sd1, sd, by = "Var2")

## this is a matrix of all distance comparisons

## need to separate out the "d vs d", "w vs w", "w vs d" - these will be the comparisons

u.sd2$Comp = paste(u.sd1$Type1,u.sd2$Type2) 

# need to replace "d w", with "w d", as these are identical
u.sd2[u.sd2 =="d w"] = "w d"

# assign the levels
u.sd2$Comp <- ordered(u.sd2$Comp,
                      levels=c("d d", "w w","w d"))
write.csv(u.sd2, here("Comp_agglom_20220411.csv"))

## use ggpubr to make a nicer boxplot with jitter datapoints

library("ggpubr")

pg = ggboxplot(u.sd2, x = "Comp", y = "value", 
          add = "jitter",
          #color = "Comp", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("d d", "w w", "w d"),
          ylab = "unweighted UniFrac distance")


#https://rpkgs.datanovia.com/ggpubr/index.html
my_comparisons = list(c("d d","w w"), c("d d", "w d"), c("w d", "w w")) # specify the comparisons

pg + stat_compare_means(comparisons = my_comparisons) # add pairwise comparisons p-value

ggsave(here("figs","20220412_comparison_unweighted_UniFrac_distance_agglomerateda.pdf"))
ggsave(here("figs","20220412_comparison_unweighted_UniFrac_distance_agglomerateda.svg"))

kruskal.test(value ~ Comp, data = u.sd2)
## there is significant difference chi-squared = 244.42, df = 2, p-value < 2.2e-16

## multiple pairwise-comparison between groups
pairwise.wilcox.test(u.sd2$value, u.sd2$Comp, p.adjust.method = "BH") #this matches the stat_compare_means

## calculate the mean and variance for each Comp group using aggregate
aggregate(u.sd2[, 3:3], list(u.sd2$Comp), mean)
aggregate(u.sd2[, 3:3], list(u.sd2$Comp), sd)

## PERMANOVA

## see for checking homogeneity condition https://mibwurrepo.github.io/Microbial-bioinformatics-introductory-course-Material-2018/beta-diversity-metrics.html
## one of the assumptions for PERMANOVA is that each group has a similar degree of multivariate scatter
sampledf = data.frame(sample_data(physeq2_red))

ps.disper = betadisper(u_u, sampledf$Type)
permutest(ps.disper,pairwise=T)
## Pr(>F) 0.001, this means the two groups do not have the same degree of multivariate scatter or multivariate spread

ps.disper = betadisper(u_u, sampledf$Location)
permutest(ps.disper,pairwise=T)

## adonis using PERMANOVA - https://github.com/joey711/phyloseq/issues/1046

adonis(u_u ~ Location, data = sampledf)
## Pr(>F) 0.371

ps.diserA = betadisper(u_u, sampledf$Tissue)
ps.diserA
#average distance to median g= 0.4553, rem = 0.4296, sg = 0.4558, w = 0.3523

permutest(ps.diserA, pairwise=T)
#Pr(>F) = 0.017 *
#g vs w = 0.002, sg vs w = 0.001, w vs. g = 0.0003, w vs. sg = 0.0004 

adonis(u_u ~ Tissue, data=sampledf)
#Tissues (df = 3), Pr(>F) = 0.001 ***

## repeat with top 500 set to see if the result is the same
u_u_500 = phyloseq::distance(physeq1_rel,m="unifrac", type=c("samples")) #unifrac is unweighted, wunifrac is weighted

u_u.m_500 = melt(as.matrix(u_u_500))

## remove self-comparisons
u_u.m_500 = u_u.m_500 %>%
  filter(as.character(Var1) != as.character(Var2)) %>%
  mutate_if(is.factor,as.character)

## order by value
u_u.m_500_2 = u_u.m_500[order(u_u.m_500$value),]

## remove every other row (these are duplicates from vice-versa comparisons) #https://statisticsglobe.com/select-odd-even-rows-columns-from-data-frame-r
row_odd <- seq_len(nrow(u_u.m_500_2)) %% 2
row_odd

u_u.m_500_3 = u_u.m_500_2[row_odd ==1, ]
dim(u_u.m_500_3) #630

sd = melt(data.frame(AC,Type))

## combined distances with sample data
## mutating joins adds coluns y to x, matching rows based on the keys (i.e., by "Var1") - https://dplyr.tidyverse.org/reference/mutate-joins.html
colnames(sd) = c("Var1", "Type1")
u.sd1_500 = left_join(u_u.m_500_3, sd, by = "Var1")

colnames(sd) = c("Var2", "Type2")
u.sd2_500 = left_join(u.sd1_500, sd, by = "Var2")

## need to separate out the "d vs d" and "w vs w" - these will be the comparisons

u.sd2_500$Comp = paste(u.sd1_500$Type1,u.sd2_500$Type2) 

# need to replace "d w", with "w d", as these are identical
u.sd2_500[u.sd2_500 =="d w"] = "w d"

u.sd2_500$Comp <- ordered(u.sd2_500$Comp,
                          levels=c("d d", "w w","w d"))

write.table(u.sd2_500,"Comptable_top500_20220414.txt", sep=";")

pgg = ggboxplot(u.sd2_500, x = "Comp", y = "value", 
          add = "jitter",
          #color = "Comp", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("d d", "w w", "w d"),
          ylab = "unweighted UniFrac distance")+
  ggtitle(paste0("top 500"))

pgg + stat_compare_means(comparisons = my_comparisons) # add pairwise comparisons p-value

ggsave(here("figs","20220414_comparison_unweighted_UniFrac_distance_top500_notaglom.pdf"))
ggsave(here("figs","20220414_comparison_unweighted_UniFrac_distance_top500_notagloma.svg"))

## calculate the mean and variance for each Comp group using aggregate
aggregate(u.sd2_500[, 3:3], list(u.sd2_500$Comp), mean)
aggregate(u.sd2_500[, 3:3], list(u.sd2_500$Comp), sd)


## kurskal-Wallis test for significance differences based on methods - non-parametric alternative to one-way ANOVA (extends Wilcoxon test)
## http://www.sthda.com/english/wiki/kruskal-wallis-test-in-r

kruskal.test(value ~ Comp, data = u.sd2_500)
## there is significant difference chi-squared = 237.55, df = 2, p-value < 2.2e-16

## multiple pairwise-comparison between groups
pairwise.wilcox.test(u.sd2_500$value, u.sd2_500$Comp,
                     p.adjust.method = "BH")

sampledf1 = data.frame(sample_data(physeq1_rel))

ps.disperB = betadisper(u_u_500, sampledf$Type)
permutest(ps.disperB,pairwise=T)

#Pr(>F) = 0.001 *** this means the two groups do not have the same degree of multivariate scatter

adonis(u_u_500 ~ Type, data = sampledf)
## Pr(>F) 0.001 ***

## results are consistent with the agglomerated results

## separate analysis for dissected and whole tissues, including identifying unique ASVs

## agglomerated merged dataset with control library removed
physeq2_red #419 ASVs

## subset only the dissected samples
dissected_split_test1 = subset_samples(physeq2_red,Type=="d")

k_test = taxa_sums(physeq2_red)
keepTaxa_test1 = taxa_names(physeq2_red)[which(k_test > 0)]
keepTaxa_test1 #409

k_test1 = taxa_sums(dissected_split_test1)
keepTaxa_dissected_test1 = taxa_names(dissected_split_test1)[which(k_test1 > 0)]
pruned_dissected_split_test1 = prune_taxa(keepTaxa_dissected_test1, dissected_split_test1) #changed from ps.filt4_merged1

## 297 taxa with 295 internal nodes and 19 samples

## ordinate using unifrac distances, weighted considers abundances, unweighted is presence/absence but no abundance)
ordu1_test1 = ordinate(pruned_dissected_split_test1, "PCoA", "unifrac")

## plot the unweighted and weighted Unifrac ordinations
diss_unifrac = plot_ordination(pruned_dissected_split_test1,ordu1_test1, type="samples", shape="Type",color="Location", title = "unweighted dissected tissues only")+
  geom_point(aes(fill=Location))+
  scale_shape_manual(values = c(16, 17,18))+
  scale_color_manual(values = c("#00AFBB","#994bc2","#ED4F50"))+
  geom_text(mapping = aes(label=Sample_AC),size=4,nudge_x=0.01,nudge_y=-0.02)+
  theme_classic() +                                                      #changes theme, removes grey background
  theme(legend.text = element_text(size = 12))+
  theme(axis.text.y.left = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.title = element_text(size=12))+
  guides(shape = guide_legend(override.aes = list(size=4)))+
  guides(fill = guide_legend(override.aes = list( size=4 )))


diss_unifrac +
  ggforce::geom_mark_ellipse(show.legend=F)+
  coord_equal()
 
ggsave(here("figs","20220414_AC_merged_unweighted_unifrac_glom_dissected.pdf"))
ggsave(here("figs","20220414_AC_merged_unweighted_unifrac_glom_dissecteda.svg"))

## similar to above, extract the distance values and then statistically compare the within and between group mean distances

## remove the rem samples (two few n=3)
no_rem = subset_samples(pruned_dissected_split_test1,Tissue!="rem")

dissected_uni = phyloseq::distance(no_rem,m="unifrac", type=c("samples")) #unifrac is unweighted, wunifrac is weighted

dissected_uni.m = melt(as.matrix(dissected_uni))

# remove self-comparisons
dissected_uni.m = dissected_uni.m %>%
  filter(as.character(Var1) != as.character(Var2)) %>%
  mutate_if(is.factor,as.character)

## order by value
dissected_uni.m2 = dissected_uni.m[order(dissected_uni.m$value),]

## remove every other row (these are duplicates from vice-versa comparisons) #https://statisticsglobe.com/select-odd-even-rows-columns-from-data-frame-r

row_odd <- seq_len(nrow(dissected_uni.m2)) %% 2
dissected_uni.m3 = dissected_uni.m2[row_odd ==1, ]
dim(dissected_uni.m3) #120

## need to add sg vs g from Type
type_ = sample_data(no_rem)$Tissue

## sample_data(no_rem)$Type = type_
AC_ = sample_data(no_rem)$Sample_AC
AC_

## need to pull from sample data vectors "Sample_AC", and "Type"
AC1 = replace(AC_, AC_=="232_sg","B5_a")
AC2 = replace(AC1, AC1=="24_g","E2_a")
AC3 = replace(AC2, AC2=="24_sg","F2_a")
AC4 = replace(AC3, AC3=="260_sg","F5_a")
AC5 = replace(AC4, AC4=="72_sg","H5_a")
AC6 = replace(AC5, AC5=="44_g","H2_c")
AC_ = AC6

sd.m = melt(data.frame(AC_,type_))

## combined distances with sample data
## mutating joins adds coluns y to x, matching rows based on the keys (i.e., by "Var1") - https://dplyr.tidyverse.org/reference/mutate-joins.html
colnames(sd.m) = c("Var1", "Type1")

dissected_uni.sd1 = left_join(dissected_uni.m3, sd.m, by = "Var1")
colnames(sd.m) = c("Var2", "Type2")
dissected_uni.sd2 = left_join(dissected_uni.sd1, sd.m, by = "Var2")

## need to separate out the "g vs g" and "sg vs g" and "sg vs sg" - these will be the comparisons

dissected_uni.sd2$Comp = paste(dissected_uni.sd1$Type1,dissected_uni.sd2$Type2) 

## 'sg g' and 'g sg' comparisons are duplicate can remove all 'g sg'
dissected_uni.sd2[dissected_uni.sd2 =="g sg"] = "sg g"

dissected_uni.sd2$Comp <- ordered(dissected_uni.sd2$Comp,
                          levels=c("g g", "sg sg","sg g"))

pgp =ggboxplot(dissected_uni.sd2, x = "Comp", y = "value", 
          add = "jitter",
          #color = "Comp", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("g g", "sg sg", "sg g"),
          ylab = "unweighted UniFrac distance")+
  ggtitle(paste0("Between and within group unweighted UniFrac distance measures - salivary and guts only"))

my_comparisons2 = list(c("g g","sg sg"), c("g g", "sg g"), c("sg sg", "sg g")) # specify the comparisons
pgp + stat_compare_means(comparisons = my_comparisons2) # add pairwise comparisons p-value
pgp + stat_compare_means(comparisons = my_comparisons2 ,method = "wilcox.test", method.args = list(alternative = "two.sided"))

###the p-values
#devtools::install_github("kassambara/ggpubr", dependencies = FALSE)
#https://github.com/kassambara/ggpubr/issues/141


ggsave(here("figs","20220414_comparison_unweighted_UniFrac_distance_guts_salivaryglandsa.pdf"))
ggsave(here("figs","20220414_comparison_unweighted_UniFrac_distance_guts_salivaryglandsb.svg"), width = 3.43, height = 4.02)

## kurskal-Wallis test for significance differences based on methods - non-parametric alternative to one-way ANOVA (extends Wilcoxon test)
## http://www.sthda.com/english/wiki/kruskal-wallis-test-in-r

kruskal.test(value ~ Comp, data = dissected_uni.sd2)
## there is slight significant difference chi-squared = 6.5856, df = 2, p-value = 0.03715

## multiple pairwise-comparison between groups
pairwise.wilcox.test(dissected_uni.sd2$value, dissected_uni.sd2$Comp,
                     p.adjust.method = "BH")

## sg vs. sg and g vs. g  p-value = 0.253, sg vs. sg and sg v g p-value = 0.470, g vs. g and sg v g p-value = 0.041

## now repeat for whole

whole_split_test1 = subset_samples(physeq2_red,Type=="w")
k_test1 = taxa_sums(whole_split_test1)
keepTaxa_whole_test1 = taxa_names(whole_split_test1)[which(k_test1 > 0)]

pruned_whole_split_test1 = prune_taxa(keepTaxa_whole_test1, whole_split_test1) 

## 257 taxa with 255 internal nodes and 17 samples

## ordinate using unifrac distances, weighted considers abundances, unweighted is presence/absence but no abundance)
ordu1_test2 = ordinate(pruned_whole_split_test1, "PCoA", "unifrac")

## plot the unweighted and weighted Unifrac ordinations
whole_unifrac = plot_ordination(pruned_whole_split_test1,ordu1_test2, type="samples",color="Location")+
  geom_point(aes(fill=Location))+
  scale_color_manual(values = c("#00AFBB","#994bc2","#ED4F50"))+
  geom_text(mapping = aes(label=Sample_AC),size=4,nudge_x=0.01,nudge_y=-0.02)+
  theme_classic() +                                                      #changes theme, removes grey background
  theme(legend.text = element_text(size = 12))+
  theme(axis.text.y.left = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.title = element_text(size=12))+
  guides(shape = guide_legend(override.aes = list(size=4)))+
  guides(fill = guide_legend(override.aes = list( size=4 )))

whole_unifrac +
  ggforce::geom_mark_ellipse(show.legend=F)+
  coord_equal()

ggsave(here("figs","20220414_Unifrac_unweighted_whole_glom.pdf"))
ggsave(here("figs","20220414_Unifrac_unweighted_whole_gloma.svg"))


##-------------------------------------------------------------------------
##                  Venn diagrams and plot unique ASVs                    -
##-------------------------------------------------------------------------

## adapted from o https://www.r-graph-gallery.com/14-venn-diagramm.html using lists of ASVs from the whole and dissected?

#install.packages("VennDiagram")
library(VennDiagram)

pruned_whole_split_test1 #257
pruned_dissected_split_test1 #297

## need to get a vectored list like before
venn_tab1 = data.frame(tax_table(pruned_dissected_split_test1)[,1:7])
dim(venn_tab1)
venn_tab2 = data.frame(tax_table(pruned_whole_split_test1)[,1:7]) #473
dim(venn_tab2)

## extract row names as vector
set_d = rownames(venn_tab1)
set_w = rownames(venn_tab2)

library(RColorBrewer)
myCol = brewer.pal(3, "Pastel2")
myCol = myCol[c(1,2)]

## Chart
venn.diagram(
  x = list(set_d, set_w),
  category.names = c("dissected" , "whole"),
  filename = '2022_02_18_venn_d_vs_wb.png',
  output=TRUE,
  
  #Output features
  #height = 480,
  #width = 480 , 
  #resolution = 300,
  #compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol, 
  
  # Numbers
  fontface = "bold",
  fontfamily = "sans"
  
)

## this venn is 112 (whole) + 145 (shared) = 257 & 152 (dissected) + 145 (shared) = 297 correct

unique_dissected1 = setdiff(set_d, set_w) #152 ASVs

only_unique_dissected1 = prune_taxa(unique_dissected1, pruned_dissected_split_test1)
only_unique_dissected1
## 152 correct

## make a table with the taxonomy
only_unique_dissected_tax_table = data.frame(tax_table(only_unique_dissected1)[,1:7])

## write this out save
write.csv(only_unique_dissected_tax_table,file=here("figs", "20220414only_unique_dissected_tax_table2.txt"),row.names = T)

## Repeat for whole
unique_whole = setdiff(set_w, set_d)
## 112 correct

## need to pull the taxonomies from the tax_tables
only_unique_whole = prune_taxa(unique_whole, pruned_whole_split_test1)
only_unique_whole_tax_table = data.frame(tax_table(only_unique_whole)[,1:7])
write.csv(only_unique_whole_tax_table,file=here("figs", "20220414only_unique_whole_tax_table.txt"),row.names = T)

## visualize at the phylum level gross differences bewteen whole and dissected unique
otu1 = otu_table(only_unique_dissected1)
taxa1 = tax_table(only_unique_dissected1)
sample1=sample_data(only_unique_dissected1)
otu2 = otu_table(only_unique_whole)
taxa2 = tax_table(only_unique_whole)
sample2=sample_data(only_unique_whole)

phy1 = phyloseq(otu_table(otu1, taxa_are_rows =FALSE),
                 sample_data(sample1),
                 tax_table(taxa1))

phy2 = phyloseq(otu_table(otu2, taxa_are_rows =FALSE),
                sample_data(sample2),
                tax_table(taxa2))

merged_unique = merge_phyloseq(phy1, phy2)
merged_unique

## 264 correct

merged_unique_prop <- transform_sample_counts(merged_unique, function(otu) otu/sum(otu))

phylumList = unique(tax_table(merged_unique_prop)[,"Phylum"])
z = dim(phylumList)

nb.cols.phy=z[1]

palette1 <- distinctColorPalette(nb.cols.phy)

## plot the relative abundance of core microbiome without controls
pp = plot_bar(merged_unique_prop, x="Sample_AC",fill = "Phylum") + 
  facet_grid(~Location1, scales="free_x",space="free_x") +
  scale_fill_manual(values = palette1)

## making a dark grey border around each ASV, with a very think line to avoid obscuring ASVs with low relative abundance
pp1 = pp + geom_bar(aes(fill=Phylum,color=""), size = .00000001,stat="identity",position = "stack") +
  scale_color_manual(values=mycolors_outline)+
  theme_bw()

## increase text size for publication
pp1 + labs(y = "Relative abundance", x = "V4 16S rDNA library")+
  theme(axis.text.x = element_text(siz=11, angle=270,vjust=0),
        axis.text.y.left = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.title = element_text(size=12),
        legend.text = element_text(size = 12)) 


ggsave(here("figs","20220414_relative_unique_phylum_merge.pdf"), width = 9, height = 6.5)
ggsave(here("figs","202200414_relative_unique_phylum_merge.svg"), width = 9, height = 6.5)

## make another venn diagram to separate based on "Tissue" (i.e., w, g, sg, w) without rem

dissected_gut = subset_samples(pruned_dissected_split_test1,Tissue=="g")
#dissected_rem = subset_samples(pruned_dissected_split_test1,Tissue=="rem")
dissected_sg = subset_samples(pruned_dissected_split_test1,Tissue=="sg")

k_test2 = taxa_sums(dissected_gut)
keepTaxa_dissected_gut1 = taxa_names(dissected_gut)[which(k_test2 > 0)]
pruned_dissected_gut1 = prune_taxa(keepTaxa_dissected_gut1, dissected_gut) #139

#k_test3 = taxa_sums(dissected_rem)
#keepTaxa_dissected_rem1 = taxa_names(dissected_rem)[which(k_test3 > 0)]
#pruned_dissected_rem1 = prune_taxa(keepTaxa_dissected_rem1, dissected_rem) #44

k_test4 = taxa_sums(dissected_sg)
keepTaxa_dissected_sg1 = taxa_names(dissected_sg)[which(k_test4 > 0)]
pruned_dissected_sg1 = prune_taxa(keepTaxa_dissected_sg1, dissected_sg) #241

pruned_whole_split_test1 #257


## need to get a vectored list like before
venn_tab3 = data.frame(tax_table(pruned_dissected_gut1)[,1:7])
#venn_tab4 = data.frame(tax_table(pruned_dissected_rem1)[,1:7])
venn_tab5 = data.frame(tax_table(pruned_dissected_sg1)[,1:7])
venn_tab2 = data.frame(tax_table(pruned_whole_split_test1)[,1:7])

## extract row names as vector
set_g = rownames(venn_tab3)
#set_rem = rownames(venn_tab4)
set_sg = rownames(venn_tab5)
set_w = rownames(venn_tab2)

library(RColorBrewer)
myCol = brewer.pal(3, "Pastel2")
myCol = myCol[c(1,2,3)]

## Chart
venn.diagram(
  x = list(set_g, set_sg, set_w),
  category.names = c("gut" , "salivary gland" , "whole"),
  filename = '2022_08_09b_venn_g_sg_w.png',
  output=TRUE,
  
  #Output features
  #height = 480,
  #width = 480 , 
  #resolution = 300,
  #compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol, 
  
  # Numbers
  fontface = "bold",
  fontfamily = "sans"
  
)


##-------------------------------------------------------------------------
##                        LDA and Lefse analysis                          -
##-------------------------------------------------------------------------

## to remove the higher taxonomic rankings in the LDA, reduce the tax_table to only a single column that includes the lowest taxonomic assignment at genus-level concatenated to the ASV identfier

## extract the taxonomy table from the ps object
tax.out = data.frame(tax_table(physeq2_red))

tax.out3 = tax.out
names = rownames(tax.out3)
tax.out3[,"ASV"] = names
tax.out3$Species_new = paste(tax.out3$Species, tax.out3$ASV, sep="_")

new_tax.out3 = dplyr::select(tax.out3, -Species, -ASV)
new_tax.out4 = dplyr::select(new_tax.out3, -Genus)
colnames(new_tax.out4)[6] <- "Genus"
new_tax.out5 = dplyr::select(new_tax.out4, -Kingdom, -Phylum, -Class, -Order, -Family)

## create a new ps object to replace with the reduced taxonomy table
physeq2_red
ps.filt2_merged_glom3 = physeq2_red
tax_table(ps.filt2_merged_glom3) <- as.matrix(new_tax.out5)

## Linear discriminant analysis (LDA) effect size (LEFSe), see: https://rdrr.io/github/yiluheihei/microbiomeMarker/f/README.md

# if (!require(remotes)) install.packages("remotes")
#remotes::install_github("yiluheihei/microbiomeMarker")
library(microbiomeMarker)

## apply the total-sum scaling 
ps.prop.filt2_merged_glom5 = transform_sample_counts(ps.filt2_merged_glom3, function(otu) otu/sum(otu))


## lefse analysis returns a microbiome biomarker stored in marker_table-class (each feature count is divied by total library size)

## Lefse analysis with microbiomeMarker -  https://yiluheihei.github.io/microbiomeMarker/articles/microbiomeMarker-vignette.html
## also see https://www.microbiomeanalyst.ca/MicrobiomeAnalyst/docs/FaqView.xhtml more explanation
## normalized total-sum scaling (TSS) followed by cumulative-sum scaling (CSS), similar approach used by Huang et al. 2021; see Paulson et al. 2013)

mm7 <- run_lefse(ps.prop.filt2_merged_glom5, norm = "CSS", group ="Type", strict = "2", sample_min = 2, lda_cutoff = 2.5)

marker_table(mm7)

## plot
plot_ef_bar(mm7, label_level=1, max_label_len = 110)+
  scale_fill_manual(values=myCol)


ggsave(here("figs","20220414_LDA_genus_by_type.pdf"), width = 7.5, height = 3.7)
ggsave(here("figs","20220414_LDA_genus_by_type.svg"), width = 7.5, height = 3.7)


mm5 <- run_lefse(ps.prop.filt2_merged_glom5, norm = "CSS", group ="Location", strict = "2", sample_min = 2, lda_cutoff = 2.5)
marker_table(mm5) # only B. miyamotoi at LP

mm6 <- run_lefse(ps.prop.filt2_merged_glom5, norm = "CSS", group ="Location1", strict = "2", sample_min = 2, lda_cutoff = 2.5)
marker_table(mm6)
myCol2 = brewer.pal(4, "Pastel2")

#plot
plot_ef_bar(mm6, label_level=1, max_label_len = 100)+
  scale_fill_manual(values=myCol2)

ggsave(here("figs","20220414_LDA_genus_by_Location1.pdf"), width = 7.5, height = 3.7)
ggsave(here("figs","20220214_LDA_genus_by_Location1.svg"), width = 7.5, height = 3.7)

#sessionInfo()

#R version 4.1.1 (2021-08-10)
#Platform: x86_64-w64-mingw32/x64 (64-bit)
#Running under: Windows 10 x64 (build 18363)

#Matrix products: default

#locale:
#  [1] LC_COLLATE=English_United States.1252 
#[2] LC_CTYPE=English_United States.1252   
#[3] LC_MONETARY=English_United States.1252
#[4] LC_NUMERIC=C                          
#[5] LC_TIME=English_United States.1252    

#attached base packages:
#  [1] grid      stats4    parallel  stats     graphics  grDevices
#[7] utils     datasets  methods   base     

#other attached packages:
#  [1] microbiomeMarker_0.99.0 RColorBrewer_1.1-2     
##[3] VennDiagram_1.6.20      futile.logger_1.4.3    
#[5] ggpubr_0.4.0            reshape2_1.4.4         
#[7] ggforce_0.3.3           phangorn_2.7.1         
#[9] ape_5.5                 DECIPHER_2.20.0        
#[11] RSQLite_2.2.8           viridis_0.6.1          
#[13] viridisLite_0.4.0       hrbrthemes_0.8.0       
#[15] scales_1.1.1            randomcoloR_1.1.0.1    
#[17] ranacapa_0.1.0          devtools_2.4.2         
#[19] usethis_2.0.1           vegan_2.5-7            
#[21] lattice_0.20-44         permute_0.9-5          
#[23] reltools_0.1.0          decontam_1.12.0        
#[25] forcats_0.5.1           stringr_1.4.0          
#[27] dplyr_1.0.7             purrr_0.3.4            
#[29] readr_2.0.1             tidyr_1.1.3            
#[31] tibble_3.1.3            ggplot2_3.3.5          
#[33] tidyverse_1.3.1         here_1.0.1             
#[35] phyloseq_1.36.0         Biostrings_2.60.2      
#[37] GenomeInfoDb_1.28.1     XVector_0.32.0         
#[39] IRanges_2.26.0          S4Vectors_0.30.0       
#[41] BiocGenerics_0.38.0    
#
#loaded via a namespace (and not attached):
#  [1] bit64_4.0.5                 knitr_1.33                 
#[3] DelayedArray_0.18.0         data.table_1.14.0          
#[5] KEGGREST_1.32.0             RCurl_1.98-1.3             
#[7] doParallel_1.0.16           generics_0.1.0             
#[9] callr_3.7.0                 microbiome_1.14.0          
#[11] lambda.r_1.2.4              bit_4.0.4                  
#[13] tzdb_0.1.2                  xml2_1.3.2                 
#[15] lubridate_1.7.10            SummarizedExperiment_1.22.0
#[17] assertthat_0.2.1            xfun_0.25                  
#[19] hms_1.1.0                   evaluate_0.14              
#[21] fansi_0.5.0                 caTools_1.18.2             
#[23] dbplyr_2.1.1                readxl_1.3.1               
#[25] igraph_1.2.6                DBI_1.1.1                  
#[27] geneplotter_1.70.0          ellipsis_0.3.2             
#[29] backports_1.2.1             V8_3.4.2                   
#[31] annotate_1.70.0             MatrixGenerics_1.4.2       
#[33] vctrs_0.3.8                 Biobase_2.52.0             
#[35] remotes_2.4.1               Cairo_1.5-12.2             
#[37] abind_1.4-5                 cachem_1.0.6               
#[39] withr_2.4.2                 treeio_1.16.2              
#[41] prettyunits_1.1.1           svglite_2.0.0              
#[43] cluster_2.1.2               lazyeval_0.2.2             
#[45] crayon_1.4.1                genefilter_1.74.0          
#[47] glmnet_4.1-2                pkgconfig_2.0.3            
#[49] labeling_0.4.2              tweenr_1.0.2               
#[51] nlme_3.1-152                pkgload_1.2.1              
#[53] rlang_0.4.11                lifecycle_1.0.0            
#[55] extrafontdb_1.0             modelr_0.1.8               
#[57] cellranger_1.1.0            rprojroot_2.0.2            
#[59] polyclip_1.10-0             matrixStats_0.60.0         
#[61] Matrix_1.3-4                aplot_0.1.0                
#[63] carData_3.0-4               Rhdf5lib_1.14.2            
#[65] reprex_2.0.1                GlobalOptions_0.1.2        
#[67] processx_3.5.2              png_0.1-7                  
#[69] rjson_0.2.20                bitops_1.0-7               
#[71] KernSmooth_2.23-20          rhdf5filters_1.4.0         
#[73] blob_1.2.2                  shape_1.4.6                
#[75] rstatix_0.7.0               gridGraphics_0.5-1         
#[77] ggsignif_0.6.2              memoise_2.0.0              
#[79] magrittr_2.0.1              plyr_1.8.6                 
#[81] gplots_3.1.1                zlibbioc_1.38.0            
#[83] compiler_4.1.1              plotROC_2.2.1              
#[85] clue_0.3-59                 DESeq2_1.32.0              
#[87] cli_3.0.1                   ade4_1.7-17                
#[89] ANCOMBC_1.2.2               patchwork_1.1.1            
#[91] ps_1.6.0                    formatR_1.11               
#[93] MASS_7.3-54                 mgcv_1.8-36                
#[95] tidyselect_1.1.1            stringi_1.7.3              
#[97] yaml_2.2.1                  locfit_1.5-9.4             
#[99] fastmatch_1.1-3             tools_4.1.1                
#[101] rio_0.5.27                  circlize_0.4.13            
#[103] rstudioapi_0.13             foreach_1.5.1              
#[105] foreign_0.8-81              gridExtra_2.3              
#[107] farver_2.1.0                Rtsne_0.15                 
#[109] digest_0.6.27               BiocManager_1.30.16        
#[111] quadprog_1.5-8              Rcpp_1.0.7                 
#[113] GenomicRanges_1.44.0        car_3.0-11                 
#[115] broom_0.7.9                 httr_1.4.2                 
#[117] gdtools_0.2.3               AnnotationDbi_1.54.1       
#[119] ComplexHeatmap_2.8.0        Wrench_1.10.0              
#[121] Rdpack_2.1.2                colorspace_2.0-2           
#[123] rvest_1.0.1                 XML_3.99-0.7               
#[125] fs_1.5.0                    splines_4.1.1              
#[127] yulab.utils_0.0.2           tidytree_0.3.5             
#[129] multtest_2.48.0             ggplotify_0.1.0            
#[131] sessioninfo_1.1.1           systemfonts_1.0.2          
#[133] xtable_1.8-4                jsonlite_1.7.2             
#[135] nloptr_1.2.2.2              ggtree_3.0.4               
#[137] futile.options_1.0.1        ggfun_0.0.3                
#[139] testthat_3.0.4              R6_2.5.1                   
#[141] pillar_1.6.2                htmltools_0.5.2            
#[143] glue_1.4.2                  fastmap_1.1.0              
#[145] BiocParallel_1.26.1         codetools_0.2-18           
#[147] pkgbuild_1.2.0              utf8_1.2.2                 
#[149] curl_4.3.2                  gtools_3.9.2               
#[151] zip_2.2.0                   openxlsx_4.2.4             
#[153] Rttf2pt1_1.3.9              limma_3.48.3               
#[155] metagenomeSeq_1.34.0        survival_3.2-12            
#[157] rmarkdown_2.10              desc_1.3.0                 
#[159] biomformat_1.20.0           munsell_0.5.0              
#[161] GetoptLong_1.0.5            rhdf5_2.36.0               
#[163] GenomeInfoDbData_1.2.6      iterators_1.0.13           
#[165] haven_2.4.3                 gtable_0.3.0               
#[167] rbibutils_2.2.3             extrafont_0.17 
