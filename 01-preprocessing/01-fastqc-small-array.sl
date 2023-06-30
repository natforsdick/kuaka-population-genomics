#!/bin/bash -e
#SBATCH -A ga03186
#SBATCH -J fastqc
#SBATCH -c 1
#SBATCH -t 00:10:00
#SBATCH -o %x.%j.%a.out
#SBATCH -e %x.%j.%a.err
#SBATCH --array=0-64%12 # Set based on filelist

cd /nesi/nobackup/ga03186/kuaka-pop-gen/output/trimmed/

SAMPLE_LIST=($(<trim-small-fastq.txt))
SAMPLE=${SAMPLE_LIST[${SLURM_ARRAY_TASK_ID}]}

module load FastQC

fastqc $SAMPLE -o ../output/trimmed-QC/

