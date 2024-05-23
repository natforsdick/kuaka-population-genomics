#!/bin/bash -e
# run as: "bash 03-filtering.sh &> filtering-out.txt"

INDIR=/nesi/nobackup/ga03186/kuaka-pop-gen/output/05-variant-calling-b/ #directory where files to filter are
vcf_out=${INDIR}filter-trial/
noLD=${vcf_out}noLD/
LD=${vcf_out}LD-filter/

STATSDIR=/nesi/nobackup/ga03186/kuaka-pop-gen/output/05-variant-calling-b/filter-trial/stats/

ml purge
ml VCFtools/0.1.15-GCC-9.2.0-Perl-5.30.1

#for loop to collect additional stats

cd $noLD

#for i in {4..5} #filtering files for 4x and 5x depth

for INBCF in *.bcf
do
    base=$(basename ${INBCF} .bcf)
    echo "Filtering SNPs for ${base}...."
    vcftools --bcf ${INBCF} \
    --out ${STATSDIR}${base} \
    --site-pi
vcftools --bcf ${INBCF} \
    --out ${STATSDIR}${base} \
    --window-pi 100000
vcftools --bcf ${INBCF} \
    --out ${STATSDIR}${base} \
    --fst-window-size 100000
vcftools --bcf ${INBCF} \
    --out ${STATSDIR}${base} \
    --SNPdensity 100000
done

cd $LD

for INBCF in *.bcf
do
    base=$(basename ${INBCF} .bcf)
    echo "Filtering SNPs for ${base}...."
    vcftools --bcf ${INBCF} \
    --out ${STATSDIR}${base} \
    --site-pi
vcftools --bcf ${INBCF} \
    --out ${STATSDIR}${base} \
    --window-pi 100000
vcftools --bcf ${INBCF} \
    --out ${STATSDIR}${base} \
    --fst-window-size 100000
vcftools --bcf ${INBCF} \
    --out ${STATSDIR}${base} \
    --SNPdensity 100000 
done
