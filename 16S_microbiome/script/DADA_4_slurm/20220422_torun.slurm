
#!/bin/bash
# SBATCH --account=def-<user> # The account to use
#SBATCH --time=08:00:00 # #The duration in HH:MM:SS format of each task in $
#SBATCH --cpus-per-task=8 #The number of cores for each task in the array
#SBATCH --mem-per-cpu=12Gb #The memory per core for each task in the array
#SBATCH --array=1           # The number of tasks:1

## set the path to the demultiplexed paired-end data
cd /global/home/hpc4605/Raw_data/Tick16sMiSeq2019/MiSeq_from_David/2021_Aug_cutadapt_demult_spacer_e0.1/Renamed/To_download_2021Aug20

## running dada2 in cac cluster

## load the r environment
module load StdEnv/2020 r/4.1.0

#R

## change jobs-script9.R to either "20220420_A_B_C_job-rscript10.R" or "20220421_A_C_job-rscript11.R"
Rscript --vanilla job-rscript9.R
