#!/bin/bash -e
#SBATCH -J merge
#SBATCH -A ga03186
#SBATCH --time=00:16:00  
#SBATCH --cpus-per-task=2
#SBATCH --array=1-69%8 # Tailor to samp number
#SBATCH --out=%x.%j.out
#SBATCH --err=%x.%j.err

###########
# PARAMS
# need to make new fq_list for BAM files
# made samplist.txt with:
# cd $BAMDIR; ls -d -1 ${PWD}/*sorted.bam > samplist.txt
BAMLIST=samplist.txt
BAMDIR=/nesi/nobackup/ga03186/kuaka-pop-gen/output/04-mapped/20241210/bam/

# MODULES
module purge
module load SAMtools/1.10-GCC-9.2.0
module list

cd $BAMDIR

QUERY1=$( awk "NR==$SLURM_ARRAY_TASK_ID" ${BAMLIST} )

echo "Task ${SLURM_ARRAY_TASK_ID}; Processing $QUERY1"

base=$(basename $QUERY1 .aligned.sorted.bam)

samtools index -b ${base}.aligned.sorted.bam

# Now let's grab some mapping stats:
echo "Getting stats"

map=$(samtools view -F4 -c ${base}.aligned.sorted.bam)
unmap=$(samtools view -f4 -c ${base}.aligned.sorted.bam)
total=$(($map + $unmap))
perc_mapped=`echo "scale=4;($map/$total)*100" | bc`

echo "$base.aligned.sorted.bam" >> ${base}-bwa_mapping_stats.txt
echo "mapped $map" >> ${base}-bwa_mapping_stats.txt
echo "% mapped $perc_mapped" >> ${base}-bwa_mapping_stats.txt
echo "unmapped $unmap" >> ${base}-bwa_mapping_stats.txt

echo "completed $QUERY1"
