#!/bin/bash -e
#SBATCH -J merge
#SBATCH -A ga03186
#SBATCH --time=00:16:00  
#SBATCH --cpus-per-task=2
#SBATCH --array=1-36%8 # Tailor to samp number
#SBATCH --out=%x.%j.out
#SBATCH --err=%x.%j.err

###########
# PARAMS
# need to make new fq_list for BAM files
BAMLIST=redo-samplist.txt
BAMDIR=/nesi/nobackup/ga03186/kuaka-pop-gen/output/04-mapped/bam/

# MODULES
module purge
module load SAMtools/1.10-GCC-9.2.0
module list

cd $BAMDIR

# first, merge sequencing duplicates from the two sequencing batches
#echo merging sequencing duplicates
#mkdir duplicates
#cp EXT041-03* duplicates/
#cp EXT041-05* duplicates/
#cp D206935_S38* duplicates/
#cp D206949_S60* duplicates/

#cd duplicates/
#samtools merge ../D206935_S38-merged.aligned.sorted.bam D206935_S38.aligned.sorted.bam EXT041-03_S3.aligned.sorted.bam
#samtools merge ../D206949_S60-merged.aligned.sorted.bam D206949_S60.aligned.sorted.bam EXT041-05_S5.aligned.sorted.bam

#cd ../
#echo completed merging

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
