#!/bin/bash -e
#SBATCH -J subsamp
#SBATCH -A ga03186
#SBATCH --time=04:30:00 # needs around 2hrs per samp pair
#SBATCH --mem=4G
#SBATCH -c 4
#SBATCH --output %x.%j.out
#SBATCH --error %x.%j.err
#SBATCH --array=0-2

# Subsampling large fastqs

ml purge
ml seqtk/1.3-gimkl-2018b

cd /nesi/nobackup/ga03186/kuaka-pop-gen/output/03-merged/

samplist=($(<subsamp-in.txt))
samp=${samplist[${SLURM_ARRAY_TASK_ID}]}

base=$(basename $samp _1.fq.gz)
echo processing ${base}

# select number of reads based on average L001_R1+L002_R1 among other samples
seqtk sample -2 -s 87 ${base}_1.fq.gz 80000000 | gzip > ${base}_sub_1.fq.gz 
seqtk sample -2 -s 87 ${base}_2.fq.gz 80000000 | gzip > ${base}_sub_2.fq.gz 

