#!/bin/bash

# SBATCH --account=def-<user>  # The account to use
#SBATCH --time=08:00:00       # The duration in HH:MM:SS format of each task in $
#SBATCH --cpus-per-task=1     # The number of cores for each task in the array
#SBATCH --mem-per-cpu=12Gb    # The memory per core for each task in the array
#SBATCH --array=1           # The number of tasks:1

cutadapt -e 0.1 --no-indels -g file:forward_barcodes.fasta -G file:reverse_barcodes.fasta -o {name1}-{name2}.1.fastq.gz -p {name1}-{name2}.2.fastq.gz sample_S1_L001_R1_001.fastq.gz sample_S1_L001_R2_001.fastq.gz