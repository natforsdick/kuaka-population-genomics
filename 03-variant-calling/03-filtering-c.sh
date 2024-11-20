#!/bin/bash -e
INDIR=/nesi/nobackup/ga03186/kuaka-pop-gen/output/05-variant-calling-b/ #directory where files to filter are
vcf_out=${INDIR}filter-trial-b/
noLD=${vcf_out}noLD-b/
LD=${vcf_out}LD-filter-b/
INBCF=Petrel_VariantCalls_concat_sort.bcf
STATSDIR=/nesi/nobackup/ga03186/kuaka-pop-gen/output/05-variant-calling-b/filter-trial-b/stats/

ml purge
ml VCFtools/0.1.15-GCC-9.2.0-Perl-5.30.1
cd $INDIR

base=$(basename ${INBCF} _concat_sort.bcf)

#for i in {4..5}
#
#do
#    vcftools --bcf ${INBCF} \
#        --out ${noLD}${base}_${i}x_coverage_0site_missing_noMinGQ.bcf \
#        --minDP ${i} \
#        --maxDP 200 \
#        --max-missing 1 \
#        --maf 0.015 \
#        --minQ 40 \
#        --remove-indels \
#        --min-alleles 2 \
#        --max-alleles 2 \
#        --remove-filtered-all \
#        --recode-bcf \
#        --recode-INFO-all
#done

ml BCFtools/1.10.2-GCC-9.2.0
for bcf in ${noLD}*0site_missing_noMinGQ*.bcf
do
    base=$(basename ${bcf} .bcf)
    echo "Running light LD pruning at 0.8 for ${base}...."
    bcftools +prune \
    -l 0.8 \
    -w 1000 \
    -O b \
    -o ${LD}${base}_0.8LD_VariantCalls.bcf \
    ${bcf}
    echo "Running moderate LD pruning at 0.6 for ${base}...."
    bcftools +prune \
    -l 0.6 \
    -w 1000 \
    -O b \
    -o ${LD}${base}_0.6LD_VariantCalls.bcf \
    ${bcf}
    echo "Running strong LD pruning at 0.4 for ${base}...."
    bcftools +prune \
    -l 0.4 \
    -w 1000 \
    -O b \
    -o ${LD}${base}_0.4LD_VariantCalls.bcf \
    ${bcf}
done

ml purge
ml VCFtools/0.1.15-GCC-9.2.0-Perl-5.30.1

for file in ${noLD}*0site_missing*.bcf
do
    base=$(basename ${file} .bcf)
    echo "Calculating depth for ${base}..."
    vcftools --bcf ${file} \
        --out ${STATSDIR}${base} \
        --site-depth
    vcftools --bcf ${file} \
        --out ${STATSDIR}${base} \
        --depth
    echo "Calculating missingness for ${base}..."
    vcftools --bcf ${file} \
        --out ${STATSDIR}${base} \
        --missing-site
    vcftools --bcf ${file} \
        --out ${STATSDIR}${base} \
        --missing-indv
    echo "Calculating individual heterozygosity for ${base}..."
    vcftools --bcf ${file} \
        --out ${STATSDIR}${base} \
        --het
done

for file in ${LD}*0site_missing*.bcf
do
    base=$(basename ${file} .bcf)
    echo "Calculating depth for ${base}..."
    vcftools --bcf ${file} \
        --out ${STATSDIR}${base} \
        --site-depth
    vcftools --bcf ${file} \
        --out ${STATSDIR}${base} \
        --depth
    echo "Calculating missingness for ${base}..."
    vcftools --bcf ${file} \
        --out ${STATSDIR}${base} \
        --missing-site
    vcftools --bcf ${file} \
        --out ${STATSDIR}${base} \
        --missing-indv
    echo "Calculating individual heterozygosity for ${base}..."
    vcftools --bcf ${file} \
        --out ${STATSDIR}${base} \
        --het
done
