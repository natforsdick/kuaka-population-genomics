#!/bin/bash -e
#SBATCH -A ga03186
#SBATCH -J pre-filt
#SBATCH -c 16
#SBATCH --mem=2G
#SBATCH --time=04:00:00 #2:30:00 for full pipeline #Walltime (HH:MM:SS)
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

# Prepping vcfs for filtering
bcfdir=/nesi/nobackup/ga03186/kuaka-pop-gen/output/05-variant-calling-b/ #bcf file directory
samplist=/nesi/nobackup/ga03186/kuaka-pop-gen/output/04-mapped/bam/samplist-2.txt

ml purge
ml BCFtools/1.15.1-GCC-11.3.0

# prepare files for filtering
for file in ${bcfdir}*.vcf
do

	base=$(basename $file .vcf)
	echo Processing ${file}
	# reheader each chunked bcf so it has the correct sample name
#	bcftools reheader -s ${samplist} ${file} -o ${bcfdir}${base}_reheader.bcf
#	wait

	#echo compressing ${file}
        bgzip -c ${bcfdir}${base}.vcf > ${bcfdir}${base}.bcf.gz
	echo indexing ${file}
	tabix ${bcfdir}${base}.bcf.gz
	# put bcf files names into a list for concatenation
	ls ${bcfdir}${base}.bcf.gz >> ${bcfdir}bcf-list.txt
done

# concatenate the chunked bcf files into a whole population bcf
echo concatenating
# -a allow overlaps, -D remove exact duplicates (outputs a single record for any duplicates
bcftools concat --file-list ${bcfdir}bcf-list.txt -a -D -O b -o ${bcfdir}Petrel_VariantCalls_concat.bcf --threads 24

