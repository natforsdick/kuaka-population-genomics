#!/bin/bash -e
#SBATCH -A ga03186
#SBATCH -J fastqc
#SBATCH -c 1
#SBATCH --mem 12G
#SBATCH -t 00:25:00
#SBATCH -o %x.%j.%a.out
#SBATCH -e %x.%j.%a.err
#SBATCH --array=0-194%10 # Set based on filelist

cd /nesi/nobackup/ga03186/kuaka-pop-gen/data/

SAMPLE_LIST=($(<mid-fastq.txt))
SAMPLE=${SAMPLE_LIST[${SLURM_ARRAY_TASK_ID}]}

module load FastQC

fastqc $SAMPLE -o ../output/raw-QC/

