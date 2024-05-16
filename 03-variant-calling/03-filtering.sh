#!/bin/bash -e
# run as: "bash 03-filtering.sh &> filtering-out.txt"

INDIR=/nesi/nobackup/ga03186/kuaka-pop-gen/output/05-variant-calling-b/ #directory where files to filter are
vcf_out=${INDIR}filter-trial/
noLD=${vcf_out}noLD/
LD=${vcf_out}LD-filter/
INBCF=Petrel_VariantCalls_concat.bcf
STATSDIR=/nesi/nobackup/ga03186/kuaka-pop-gen/output/05-variant-calling-b/filter-trial/stats/

ml purge 
ml VCFtools/0.1.15-GCC-9.2.0-Perl-5.30.1

cd $INDIR
mkdir -p $vcf_out $noLD $LD $STATSDIR

#for loop to filter file with different values for parameters including
#missingness, depth, and GQ

base=$(basename ${INBCF} _concat.bcf)

#for i in {4..5} #filtering files for 4x and 5x depth

for i in 5
do
    echo "Filtering SNPs for ${base}...." 
    vcftools --bcf ${INBCF} \
        --out ${noLD}${base}_${i}x_coverage_0site_missing_noMinGQ.bcf \
        --minDP ${i} \
        --maxDP 200 \
        --max-missing 1 \
        --maf 0.05 \
        --minQ 20 \
        --remove-indels \
	--min-alleles 2 \
	--max-alleles 2 \
        --remove-filtered-all \
        --recode-bcf \
        --recode-INFO-all 
    vcftools --bcf ${INBCF} \
        --out ${noLD}${base}_${i}x_coverage_0.1site_missing_noMinGQ.bcf \
        --minDP ${i} \
        --maxDP 200 \
        --max-missing 0.9 \
        --maf 0.05 \
        --minQ 20 \
        --min-alleles 2 \
        --max-alleles 2 \
        --remove-indels \
        --remove-filtered-all \
        --recode-bcf \
        --recode-INFO-all 
    vcftools --bcf ${INBCF} \
        --out ${noLD}${base}_${i}x_coverage_0.2site_missing_noMinGQ.bcf \
        --minDP ${i} \
        --maxDP 200 \
        --max-missing 0.8 \
        --maf 0.05 \
        --minQ 20 \
        --min-alleles 2 \
        --max-alleles 2 \
        --remove-indels \
        --remove-filtered-all \
        --recode-bcf \
        --recode-INFO-all 
    vcftools --bcf ${INBCF} \
        --out ${noLD}${base}_${i}x_coverage_0site_missing_MinGQ10.bcf \
        --minDP ${i} \
        --maxDP 200 \
        --max-missing 1 \
        --maf 0.05 \
        --minQ 20 \
        --minGQ 10 \
        --min-alleles 2 \
        --max-alleles 2 \
        --remove-indels \
        --remove-filtered-all \
        --recode-bcf \
        --recode-INFO-all 
    vcftools --bcf ${INBCF} \
        --out ${noLD}${base}_${i}x_coverage_0.1site_missing_MinGQ10.bcf \
        --minDP ${i} \
        --maxDP 200 \
        --max-missing 0.9 \
        --maf 0.05 \
        --minQ 20 \
        --minGQ 10 \
        --min-alleles 2 \
        --max-alleles 2 \
        --remove-indels \
        --remove-filtered-all \
        --recode-bcf \
        --recode-INFO-all 
    vcftools --bcf ${INBCF} \
        --out ${noLD}${base}_${i}x_coverage_0.2site_missing_MinGQ10.bcf \
        --minDP ${i} \
        --maxDP 200 \
        --max-missing 0.8 \
        --maf 0.05 \
        --minQ 20 \
        --minGQ 10 \
        --min-alleles 2 \
        --max-alleles 2 \
        --remove-indels \
        --remove-filtered-all \
        --recode-bcf \
        --recode-INFO-all 
    vcftools --bcf ${INBCF} \
        --out ${noLD}${base}_${i}x_coverage_0site_missing_MinGQ20.bcf \
        --minDP ${i} \
        --maxDP 200 \
        --max-missing 1 \
        --maf 0.05 \
        --minQ 20 \
        --minGQ 20 \
        --min-alleles 2 \
        --max-alleles 2 \
        --remove-indels \
        --remove-filtered-all \
        --recode-bcf \
        --recode-INFO-all 
    vcftools --bcf ${INBCF} \
        --out ${noLD}${base}_${i}x_coverage_0.1site_missing_MinGQ20.bcf \
        --minDP ${i} \
        --maxDP 200 \
        --max-missing 0.9 \
        --maf 0.05 \
        --minQ 20 \
        --minGQ 20 \
        --min-alleles 2 \
        --max-alleles 2 \
        --remove-indels \
        --remove-filtered-all \
        --recode-bcf \
        --recode-INFO-all 
    vcftools --bcf ${INBCF} \
        --out ${noLD}${base}_${i}x_coverage_0.2site_missing_MinGQ20.bcf \
        --minDP ${i} \
        --maxDP 200 \
        --max-missing 0.8 \
        --maf 0.05 \
        --minQ 20 \
        --minGQ 20 \
        --min-alleles 2 \
        --max-alleles 2 \
        --remove-indels \
        --remove-filtered-all \
        --recode-bcf \
        --recode-INFO-all
done
wait

ml BCFtools/1.10.2-GCC-9.2.0

echo "Filtering for Linkage parameters..."
#for loop to filter previous filtered files for linkage
for bcf in ${noLD}*.bcf
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
wait

ml purge
ml VCFtools/0.1.15-GCC-9.2.0-Perl-5.30.1

#calculating statistics for no linkage filtered files
for file in ${noLD}*.bcf
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
wait

#calculating statistics for linkage filtered files
for file in ${LD}*.bcf
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
wait
