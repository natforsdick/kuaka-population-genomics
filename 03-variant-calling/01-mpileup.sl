#!/bin/bash -e
#SBATCH -A ga03186
#SBATCH -J mpileup
#SBATCH -c 32
#SBATCH --mem=12G
#SBATCH --time=06:00:00 #Walltime (HH:MM:SS) 
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

ref=/nesi/nobackup/ga03186/kuaka-genome/05-scaffolding/05_yahs/01-kuaka-hifiasm-p_ctg-purged-yahs_scaffolds_final.fa
bamdir=/nesi/nobackup/ga03186/kuaka-pop-gen/output/04-mapped/bam/
bcfdir=/nesi/nobackup/ga03186/kuaka-pop-gen/output/05-variant-calling/ # output bcf file directory
samplist=/nesi/nobackup/ga03186/kuaka-pop-gen/output/04-mapped/bam/samplist.txt
platform="Illumina"

ml purge
ml BCFtools/1.15.1-GCC-11.3.0 SAMtools/1.15.1-GCC-11.3.0

cd $bamdir 

if [ ! -e ${samplist} ]; then
	ls -d ${bamdir}*.bam > samplist.txt
fi

if [ ! -e ${bcfdir}chunks/ ]; then
	mkdir -p ${bcfdir}chunks
fi

#chunk bam files into 12 pieces using custom perl script 
#@Lanilen/SubSampler_SNPcaller/split_bamfiles_tasks.pl
#chunked files will help mpileup run faster
echo chunking
perl /nesi/project/ga03186/kuaka-pop-gen/03-variant-calling/split_bamfiles_tasks.pl -b ${samplist} -g $ref -n 12 -o ${bcfdir}chunks | parallel -j 12 {}

#echo chunking complete
#run bcftools mpileup in parallel on chunks of bam files with BCFtools v 1.11
echo running mpileup
for (( i=1; i<=12; i++ )); do
        bcftools mpileup -E -O b -f $ref -a AD,ADF,DP,ADR,SP -o ${bcfdir}kuaka_${i}_raw.bcf ${bcfdir}chunks/${i}/* &
done
wait
echo mpileup complete

#SNP calling on bcf files with bcftools call
for file in ${bcfdir}*.bcf
	do
	base=$(basename $file .bcf)
	bcftools call $file --threads 24 -mv -O v -o ${bcfdir}${base}_VariantCalls.vcf &
done
wait
echo variant calling is complete
