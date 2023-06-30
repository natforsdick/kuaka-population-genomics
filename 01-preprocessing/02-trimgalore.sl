#!/bin/bash -e
#SBATCH -A ga03186
#SBATCH -J trimgalore
#SBATCH -c 4
#SBATCH --mem=1G
#SBATCH --array=0
#SBATCH --time=02:30:00 #Walltime (HH:MM:SS) 
#SBATCH --output=%x.%j.%a.out
#SBATCH --output=%x.%j.%a.err

# PARAMS
INDIR=/nesi/nobackup/ga03186/kuaka-pop-gen/data/
OUTDIR=/nesi/nobackup/ga03186/kuaka-pop-gen/output/02-trimmed/02b-trimgalore/
SAMPLE_LIST=($(<${INDIR}mid-fastq1.txt))
SAMPLE=${SAMPLE_LIST[${SLURM_ARRAY_TASK_ID}]}

base=$(basename $SAMPLE _R1_001.fastq.gz)

#trim files with TrimGalore v. 0.6.6 before alignment- change parameters based on results
ml purge
ml TrimGalore/0.6.7-gimkl-2020a-Python-3.8.2-Perl-5.30.1

cd $INDIR

trim_galore --paired --nextseq 28 --length 50 \
--three_prime_clip_R1 5 --three_prime_clip_R2 5 --clip_R1 20 --clip_R2 20 --2colour 20 \
--fastqc ${INDIR}${base}_R1_001.fastq.gz ${INDIR}${base}_R2_001.fastq.gz \
-o ${OUTDIR} --basename ${base}
# trims paired end reads of a sample to a minimum length of 50, does a 3' clip of 5 bp and 5' clip of 20 bp,
# performs clips with a 2-colour compatible quality phred score of 20, and clip of 20 bases

wait
echo "done trimming ${base}"

