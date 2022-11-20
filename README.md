---
title: "README.Rmd"
author: "Amber Paulson"
date: "04/01/2022"
output: html_document
runtime: shiny
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

<br />
<div align="center">
  <h3 align="center">Tick_microbiome README </h3>
  <p align="center">
    <a href="https://github.com/damselflywingz/tick_microbiome/">github.com/damselflywingz/tick_microbiome</a>
    <br />
  </p>
</div>



<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li><a>[Introduction](#introduction)</a></li>
    <li><a>[Built With](#built-with)</a></li>
    <li><a>[Getting Started](#getting-started)</a></li>
    <li><a>[Demultiplex](#demultiplex)</a></li>
    <li><a>[Filter Denoise](#filter-denoise)</a></li>
    <li><a>[Decontaminate](#decontaminate)</a></li>
    <li><a>[Statistical Analysis](#statistical-analysis)</a></li>
    <li><a>[BLAST/ Phylogeny](#blast-phylogeny)</a></li>
    <li><a>[Contact](#contact)</a></li>
    </ol>
</details>

### Introduction

This README provides scripts and other resources used to analyze 16S rRNA sequencing from the microbiome of the blacklegged tick, <i>Ixodes scapularis</i> - see Preprint Paulson A.R., Lougheed, S.C., Huang, D., and Colautti, R.I. 2022. Multi-omics analysis identifies symbionts and pathogens of blacklegged ticks (<i>Ixodes scapularis</i>) from a Lyme disease hotspot in southeastern Ontario, Canada.

https://doi.org/10.1101/2022.11.09.515820

<p align="right">(<a href="#top">back to top</a>)</p>

### Built With 

* <a href="https://cutadapt.readthedocs.io/en/stable/?target=_blank" target="_blank">Cutadapt</a>
* <a href="https://benjjneb.github.io/dada2/index.html" target="_blank">DADA2</a>
* <a href="https://benjjneb.github.io/decontam/" target="_blank">decontam</a>
* <a href="https://joey711.github.io/phyloseq/" target="_blank">phyloseq</a>

<p align="right">(<a href="#top">back to top</a>)</p>


<!-- GETTING STARTED -->
## Getting Started

The analysis involves:
<ol>
    <li>Demultiplexing raw data using Cutadapt Version 2.6;</li>
<li>Filtering, trimming, denoising, merging, collapsing, chimera removal, and taxonomic assignment using DADA2 Version 1.20.0;</li>
<li>Contamination removal using decontam Version 1.12.0;</li>
<li>Statistical analysis using phyloseq Version 1.36.0;</li>
<li>local BLAST searches using rBLAST Version 0.99.2; and</li>
<li>Phylogenetic reconstruction using phangorn Version 2.7.1.</li>
</ol>

The first two steps of this pipeline (e.g., demultiplexing, taxonomic assignment, denoising, etc.) involve computationally intensive processes, so these were computed on the Queen's University Center for Advanced Computing's Frontenac high-performance computing cluster. 

Since the last two steps of this pipeline (i.e., decontaminating and statistical analysis) are less computationally intensive, these steps were executed locally using R in RStudio (noting phylogeny building with 500 sequences on a single core can take about ~ 30 min).

<p align="right">(<a href="#top">back to top</a>)</p>



<!-- DEMULTIPLEX -->
## Demultiplex

Starting with the raw non-demultiplexed sequence data [SRR17194087](https://www.ncbi.nlm.nih.gov/sra/SRR17194087), the non-interleaved raw reads must be sorted with barcodes and spacers in combinatorial format using Cutadapt. All the required files are located in the [github.com/damselflywingz/tick_microbiome/tree/main/cutadapt/](https://github.com/damselflywingz/tick_microbiome/tree/main/cutadapt/) directory, please see there for more information.

To install Cutadapt follow the instructions to create a Conda environment, and then activate the cutadapt environment - [cutadapt.readthedocs.io/en/stable/installation](https://cutadapt.readthedocs.io/en/stable/installation.html).

The task was submitted through the SLURM scheduler (see: [cutadapt_demultiplex.slurm](https://github.com/damselflywingz/tick_microbiome/tree/main/cutadapt/cutadapt_demultiplex.slurm)).

After the job is complete the mmv command is used to batch rename the files

<p style="font-family: Lucida Console, monospace; font-size:11pt;background-color: #DDDDDD">
mmv < patterns_update.txt
</p>

<p align="right">(<a href="#top">back to top</a>)</p>


<!-- FILTER DENOISE -->
## Filter Denoise

Before proceeding please see [github.com/damselflywingz/tick_microbiome/tree/main/16S_microbiome/script/DADA_4_slurm/Not_sequenced_note.txt](https://github.com/damselflywingz/tick_microbiome/blob/main/16S_microbiome/script/DADA_4_slurm/Not_sequenced_note.txt) as some of the libraries were not sequenced and need to be removed prior to running the pipeline. These files need to be moved into a "Not_sequenced folder".

Working with demultiplexed paired-end libraries, the DADA2 processing follows the tutorial at <a href="https://benjjneb.github.io/dada2/index.html" target="_blank">benjjneb.github.io/dada2/tutorial </a> for the most part. See the link for more information and resources.


There are two separate R scripts:
<ol>
    <li>[20220420_A_B_C_job-rscript10.R](https://github.com/damselflywingz/tick_microbiome/blob/main/16S_microbiome/script/DADA_4_slurm/20220420_A_B_C_job-rscript10.R); and </li>
    <li>[20220421_A_C_job-rscript11.R](https://github.com/damselflywingz/tick_microbiome/blob/main/16S_microbiome/script/DADA_4_slurm/20220421_A_C_job-rscript11.R).</li>
    </ol>

The script requires the RDP training set and species-specific reference databases to be available in the "reference" directory, see: [github.com/damselflywingz/tick_microbiome/tree/main/16S_microbiome/reference](https://github.com/damselflywingz/tick_microbiome/tree/main/16S_microbiome/reference). There are other reference test sets available <a href="https://benjjneb.github.io/dada2/training.html" target="_blank">benjjneb.github.io/dada2/training </a>for more options and the latest versions of the training sets. 

The DADA2 R scripts can be submitted using  [20220422_torun.slurm](https://github.com/damselflywingz/tick_microbiome/blob/main/16S_microbiome/script/DADA_4_slurm/20220422_torun.slurm)

<p align="right">(<a href="#top">back to top</a>)</p>
 

<!-- DECONTAMINATE -->
## Decontaminate

The decontamination step (and the following statistical analysis step) were completed locally with RStudio Version 1.2.5033 and R version 4.1.0.

First download the entire [github.com/damselflywingz/tick_microbiome/tree/main/16S_microbiome/](https://github.com/damselflywingz/tick_microbiome/tree/main/16S_microbiome/) directory and contents.

Then load the [16S_microbiome.Rproj](https://github.com/damselflywingz/tick_microbiome/tree/main/16S_microbiome/16S_microbiome.Rproj) in RStudio. 

Then open the R script [github.com/damselflywingz/tick_microbiome/tree/main/16S_microbiome/script/16S_decontam_analysis_Paulson_v20220401.R](https://github.com/damselflywingz/tick_microbiome/tree/main/16S_microbiome/script/16S_decontam_analysis_Paulson_v20220401.R). 

<p align="right">(<a href="#top">back to top</a>)</p>



<!-- STATISTCAL ANALYSIS -->
## Statistical Analysis

Following the decontaminate step (see above), the statistical analysis was completed by continuing through the [16S_decontam_analysis_Paulson_v20220401.R](https://github.com/damselflywingz/tick_microbiome/tree/main/16S_microbiome/script/16S_decontam_analysis_Paulson_v20220401.R) script. See the script for more detailed information.

<p align="right">(<a href="#top">back to top</a>)</p>



<!-- BLAST/ PHYLOGENY -->
## BLAST/ Phylogeny

Local BLAST search of the GenBank's 16S reference database and phylogenetic analysis of selected V4 16S rDNA sequences from the tick microbiome were completed by following the BlastPhylo script created by David Huang. [BlastPhylo_v20220418.R](https://github.com/damselflywingz/tick_microbiome/blob/main/16S_microbiome/script/BlastPhylo_v20220418.R). See the script for more detailed information.

<p align="right">(<a href="#top">back to top</a>)</p>



<!-- CONTACT -->
## Contact

Amber Paulson - <a href="https://twitter.com/dragonflywingz" target="_blank">@Dragonflywingz </a> - Amber[dot]Rose[dot]Paulson[at]gmail[dot]com

Project Link: [https://github.com/damselflywingz/tick_microbiome](https://github.com/damselflywingz/tick_microbiome)

Licensed under Creative Commons Attribution 4.0 International (CC-BY 4.0) [https://creativecommons.org/licenses/by/4.0/](https://creativecommons.org/licenses/by/4.0/)

<p align="right">(<a href="#top">back to top</a>)</p>

