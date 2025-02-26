#!/bin/bash -e
#SBATCH -A ga03186
#SBATCH -J pre-filt
#SBATCH -c 18
#SBATCH --mem=8G
#SBATCH --time=10:00:00 #2:30:00 for full pipeline #Walltime (HH:MM:SS)
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

# Prepping vcfs for filtering
bcfdir=/nesi/nobackup/ga03186/kuaka-pop-gen/output/05-variant-calling-20241210/ #bcf file directory

ml purge
ml BCFtools/1.15.1-GCC-11.3.0

# prepare files for filtering
for file in ${bcfdir}*VariantCalls.vcf
do

	base=$(basename $file .vcf)
	echo Processing ${file}

        bgzip -c ${bcfdir}${base}.vcf > ${bcfdir}${base}.bcf.gz
	echo indexing ${file}
	tabix ${bcfdir}${base}.bcf.gz
	# put bcf files names into a list for concatenation
	ls ${bcfdir}${base}.bcf.gz >> ${bcfdir}bcflist.txt
done

# concatenate the chunked bcf files into a whole population bcf
echo concatenating
# -a allow overlaps, -D remove exact duplicates (outputs a single record for any duplicates)
bcftools concat --file-list ${bcfdir}bcflist.txt -a -D -O b -o ${bcfdir}Petrel_VariantCalls_concat.bcf --threads 32

# now let's sort the resultant bcf
ml purge; ml BCFtools/1.19-GCC-11.3.0
TMPDIR=${bcfdir}temp
mkdir -p $TMPDIR

bcftools sort --temp-dir $TMPDIR -O b -o ${bcfdir}Petrel_VariantCalls_concat_sort.bcf ${bcfdir}Petrel_VariantCalls_concat.bcf
