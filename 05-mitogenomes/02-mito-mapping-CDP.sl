#!/bin/bash -e
#SBATCH -A ga03186
#SBATCH -J mitogenome-mapping
#SBATCH --cpus-per-task=6
#SBATCH --mem=22G 
#SBATCH -t 4:00:00 # initially 2:20 to get through the bulk of the samples (excluding decompression) - up to 8 hrs for large input files
#SBATCH --out %x.%j.%a.out
#SBATCH --err %x.%j.%a.err
#SBATCH --array=1-27%8 # test with a small handful to start

# Extracting mitogenomes from short-read WGS data
# This script expects a make_coverage_plots Rscript to be in the processing directory
# This script processes multiple paired-end fastq files through to merged bam files

####################
# variables
# reference file
refdir=/nesi/nobackup/ga03186/kuaka-genome/mitohifi/
reffile=NC_052809.1-P-urinatrix-mitogenome-50bprem.fasta
ref=$refdir$reffile

# path to trimmed short-read data
indir=/nesi/nobackup/ga03186/kuaka-pop-gen/output/03-merged/
outdir=/nesi/nobackup/ga03186/kuaka-pop-gen/output/mitogenomes/CDP/
cd $outdir

# samplist containing list of file prefixes for the trimmed short-read data
samplist=${indir}mitolist-CDP.txt
fq1=_1.fq
fq2=_2.fq

# the name of the target mitosequence 
mtchr="NC_052809.1"
platform="Illumina"
library="Lib1" 

####################
# software
module purge
module load picard/2.26.10-Java-11.0.4 BWA/0.7.18-GCC-12.3.0 

####################
# prep to run as array
samp=$( awk "NR==$SLURM_ARRAY_TASK_ID" ${samplist} )
echo "Task ${SLURM_ARRAY_TASK_ID}; Processing $samp"

####################
echo "processing data for $samp"
# unzip data if necessary
if [ -e $indir${samp}${fq1}.gz ] || [ -e $indir${samp}${fq2}.gz ]; then
    echo "Unzipping fastq files"
    gunzip $indir${samp}${fq1}.gz
    gunzip $indir${samp}${fq2}.gz
fi

#get readgroup info
##split the machine/lane info from the fq file
infoline=$(head -n 1 ${indir}$samp$fq1)
instrument=`echo $infoline | cut -d ':' -f1`
instrument=$(echo $instrument | sed -e 's/ /_/g') # this is to fix issues with spaces in the instrument name
instrumentrun=`echo $infoline | cut -d ':' -f2` 
flowcell=`echo $infoline | cut -d ':' -f3`
lane=`echo $infoline | cut -d ':' -f4`
index=`echo $infoline | cut -d ':' -f10`
#work that info into some 
rgid="ID:${instrument}_${instrumentrun}_${flowcell}_${lane}_${index}"
rgpl="PL:${platform}"
rgpu="PU:${flowcell}.${lane}"
rglb="LB:${samp}_${library}"
rgsm="SM:${samp}"
echo $rgid

#specify the paired-end reads
alignfastq1=${indir}$samp$fq1
alignfastq2=${indir}$samp$fq2

# align paired-end files with bwa mem
echo "aligning and sorting data for $samp"
bwa mem -t 16 -M -R @RG'\t'$rgid'\t'$rgpl'\t'$rgpu'\t'$rgsm $ref $alignfastq1 $alignfastq2 | picard SortSam -I /dev/stdin -O ${outdir}${samp}_sort.bam --SORT_ORDER coordinate --TMP_DIR ./tmp
module purge; module load SAMtools/1.19-GCC-12.3.0 

# index the sorted BAM-formatted alignments
echo "indexing the bamfile for $samp"
alignedbam=${outdir}${samp}_sort.bam
samtools index -@ 12 -b $alignedbam

# get stats from the indexed BAM file
echo "collecting stats for $samp"
samtools stats -@ 12 $alignedbam >  ${outdir}stats/${samp}_sort_stats.txt
samtools idxstats -@ 12 $alignedbam > ${outdir}stats/${samp}_sort_idxstats.txt

# clean the bam file: Fixes errant MAPQ scores (unmapped reads with scores not zero), soft clips pairs hanging off alignments
echo "cleaning the bamfile $samp"
module purge; module load picard/2.26.10-Java-11.0.4
picard CleanSam -I $alignedbam -O ${outdir}${samp}_allsortclean.bam --CREATE_INDEX true

# remove unmapped reads
echo "removing unmapped reads for $samp"
module purge; module load SAMtools/1.19-GCC-12.3.0 
samtools view -@ 12 -b -F 0x0004 ${outdir}${samp}_allsortclean.bam > ${outdir}${samp}_maponly.bam

# remove duplicates 
echo "removing duplicates for $samp"
module purge; module load picard/2.26.10-Java-11.0.4
predup=${outdir}${samp}_maponly.bam
picard MarkDuplicates \
    -I $predup \
    -O ${outdir}${samp}_remdup.bam \
    --METRICS_FILE ${outdir}stats/${samp}_remdup_metrics.txt \
    --CREATE_INDEX true \
    --ASSUME_SORTED true
echo "finished processing $samp"

####################
echo "making plots and stats for $samp"
#create coverage plot
if [ -e /nesi/project/ga03186/kuaka-pop-gen/05-mitogenomes/make_coverage_plots2.0.R ]; then
	echo "Creating coverage plot"
	bamprefix="_remdup.bam"
    cd ${outdir}coverage/
	# create coverage file
    module purge; module load SAMtools/1.19-GCC-12.3.0
	samtools depth -r $mtchr -a ${outdir}${samp}${bamprefix} > ${outdir}coverage/${samp}_coverage.txt
	# create coverage plot
    module purge; module load R/4.3.1-gimkl-2022a
	Rscript /nesi/project/ga03186/kuaka-pop-gen/05-mitogenomes/make_coverage_plots2.0.R $samp $mtchr
	cd ..
else
    echo "Cannot find make_coverage_plots.R: skipping coverage plot step"
fi
cd ..
echo "completed processing $samp"
####################
