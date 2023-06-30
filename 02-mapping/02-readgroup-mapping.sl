#!/bin/bash -e
#SBATCH -J bwa
#SBATCH -A ga03186
#SBATCH --time=02:30:00 
#SBATCH --mem=28G 
#SBATCH --cpus-per-task=32
#SBATCH --array=1-2 #14%7 # Tailor to samp number
#SBATCH --out=%x.%j.out
#SBATCH --err=%x.%j.err


# MODULES
module purge
module load  BWA/0.7.17-GCC-9.2.0 SAMtools/1.10-GCC-9.2.0
module list

###########

# PARAMS
# need to make new fq_list for trimgalore outputs
fq_list=/nesi/nobackup/ga03186/kuaka-pop-gen/output/03-merged/merged-samplist.txt

reffile=01-kuaka-hifiasm-p_ctg-purged-yahs_scaffolds_final
refdir=/nesi/nobackup/ga03186/kuaka-genome/05-scaffolding/05_yahs/
ref=$refdir$reffile

INDIR=/nesi/nobackup/ga03186/kuaka-pop-gen/output/03-merged/
SAMDIR=/nesi/nobackup/ga03186/kuaka-pop-gen/output/04-mapped/sam/
BAMDIR=/nesi/nobackup/ga03186/kuaka-pop-gen/output/04-mapped/bam/
fq1=_val_1.fq.gz #Read 1 suffix
fq2=_val_2.fq.gz #Read 2 suffix
# fq for subsampled large files
#fq1=_val_sub_1.fq.gz
#fq2=_val_sub_2.fq.gz

platform="Illumina"

###########

# MAPPING
cd $INDIR

QUERY1=$( awk "NR==$SLURM_ARRAY_TASK_ID" ${fq_list} ) #`cat ${fq_list} | awk -v line=$SLURM_ARRAY_TASK_ID '{if(NR == line) print $1}'`

echo "Task ${SLURM_ARRAY_TASK_ID}; Processing $QUERY1"

# capture readgroup info
base=$(basename $QUERY1 _val_1.fq.gz)
infoline=$(zcat ${QUERY1} | head -n 1)
instrument=`echo ${infoline} | cut -d ':' -f1`
instrumentrun=`echo $infoline | cut -d ':' -f2`
flowcell=`echo $infoline | cut -d ':' -f3`
lane=`echo $infoline | cut -d ':' -f4`
index=`echo $infoline | cut -d ':' -f10`

# incorporate sample information into the alignment
rgid="ID:${instrument}_${instrumentrun}_${flowcell}_${lane}_${index}"
rgpl="PL:${platform}"
rgpu="PU:${flowcell}.${lane}"
rglb="LB:${base}_library1"
rgsm="SM:${base}"

echo "Aligning reads for $base"
bwa mem -M -R @RG'\t'$rgid'\t'$rgpl'\t'$rgpu'\t'$rglb'\t'$rgsm -t 32 $ref $QUERY1 ${datadir}${base}${fq2} > ${SAMDIR}${base}.sam

echo "Converting sam file to bam file for $base"
samtools view -@ 32 -T $ref.fa -b ${SAMDIR}${base}.sam > ${BAMDIR}${base}.bam

echo "Sorting and indexing file"
samtools sort -@ 32 -o ${BAMDIR}${base}.aligned.sorted.bam ${BAMDIR}${base}.bam

echo "Removing intermediate files"
rm ${SAMDIR}${base}.sam
rm ${BAMDIR}${base}.bam
echo "Completed $QUERY1"
