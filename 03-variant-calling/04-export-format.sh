#!/bin/bash -e

cd /nesi/nobackup/ga03186/kuaka-pop-gen/output/05-variant-calling/filter-trial/
mkdir -p ./export/

ml purge
ml BCFtools/1.10.2-GCC-9.2.0 Stacks/2.61-gimkl-2022a

bcftools convert -O v -o ./noLD/Petrel_VariantCalls_4x_coverage_0.2site_missing_MinGQ10.bcf.recode.vcf ./noLD/Petrel_VariantCalls_4x_coverage_0.2site_missing_MinGQ10.bcf.recode.bcf
bcftools convert -O v -o ./LD-filter/Petrel_VariantCalls_5x_coverage_0.1site_missing_MinGQ10.bcf.recode_0.4LD_VariantCalls.vcf LD-filter/Petrel_VariantCalls_5x_coverage_0.1site_missing_MinGQ10.bcf.recode_0.4LD_VariantCalls.bcf
populations -V ./noLD/Petrel_VariantCalls_4x_coverage_0.2site_missing_MinGQ10.bcf.recode.vcf -O ./export/ -M popmap-2.txt --plink
populations -V ./LD-filter/Petrel_VariantCalls_5x_coverage_0.1site_missing_MinGQ10.bcf.recode_0.4LD_VariantCalls.vcf -O ./export/ -M popmap-2.txt --plink

ml purge
ml PLINK/1.09b6.16

plink --file ./export/Petrel_VariantCalls_4x_coverage_0.2site_missing_MinGQ10.bcf.recode.p.plink --aec --recode A
plink --file ./export/Petrel_VariantCalls_5x_coverage_0.1site_missing_MinGQ10.bcf.recode_0.4LD_VariantCalls.p.plink --aec --recode A
