#!/bin/bash -e
#SBATCH -A ga03186
#SBATCH -J vcfcall
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH -t 1:00:00
#SBATCH --out %x.%j.out
#SBATCH --err %x.%j.err

#call variants from BAM files produced by 02-mito-mapping.sl
#Uses GATK's HaploCaller to call all variants, and combine into a group VCF 
# Protocol by Sophia Cameron-Christie, modified by Nat Forsdick

####################
#variables 
# the name and location of the target mitosequence 
refdir=/nesi/nobackup/ga03186/kuaka-genome/mitohifi/
reffile=kuaka_final_mitogenome-80bprem.fasta
ref=$refdir$reffile

# path to trimmed shot-read data
indir=/nesi/nobackup/ga03186/kuaka-pop-gen/output/03-merged/
outdir=/nesi/nobackup/ga03186/kuaka-pop-gen/output/mitogenomes/KWH/

grpname=KWH-mito #name of this group of samples
projectdir=/nesi/nobackup/ga03186/kuaka-pop-gen/output/mitogenomes/
grpcalldir="call_group_genotypes"

ploidy=1 #the ploidy of the data (eg, 1 for mitochondria, 2 for autosomal)
callingtype="split"
filtertype="hard"
# How to filter the variants; either "VQSR" or "hard".
# For genome-wide data, the recommendation is Variant Quality Score Recalibration (VQSR) but this is not recommended for when only mitochondrial data is available
# because the number of variants are too few; therefore "hard" filtering is used for these variants.
# Hard filters SHOULD BE EXAMINED AND ADAPTED according to the dataset; further discussion can be found on GATK's site: software.broadinstitute.org/gatk/guide/article?id=6925

MQthreshold="40.0" #the quality threshold for filtering - should be specific to the group # set to 40 initially, can reduce this if required.
call="yes" #either yes or [anything else]. If "yes", variants will be called from sample BAMs. If it is not yes, existing sample VCFs will be used.

# a list containing sample names as they were processed in 02-mito-mapping.sl
samplist=$(cat ${indir}mitolist-KWH.txt)
vcflist=${outdir}KWH2-vcf.list
####################

cd $outdir

#specify the name of the BAM file to call from
bamprefix="_remdup.bam"

module purge
module load SAMtools/1.19-GCC-12.3.0 picard/2.26.10-Java-11.0.4

if [ ! -d "$grpcalldir" ]; then
echo "creating new directory for group variant calling"
mkdir $grpcalldir
fi

idxfile=${reffile%.*}.fai
if [ ! -e "$refdir$reffile.fai" ]; then
echo "Reference is not indexed: indexing with SAMTOOLS"
samtools faidx $ref
else
echo "Index for reference file found"
fi

dictfile=${reffile%.*}.dict
if [ ! -e "$refdir${reffile%.*}.dict" ]; then
echo "Dictionary file of reference does not exist: creating dictionary"
picard CreateSequenceDictionary -R $ref -O $refdir${reffile%.*}.dict
else
echo "Dictionary file found"
fi

cd $grpcalldir

if [ $callingtype = "split" ]; then

#begin the calling on each sample
if [ $call = "yes" ]; then
echo "variants will be called from individual sample BAMs"
for samp in $samplist
do

# call variants for sample
module purge; module load GATK/4.5.0.0-gimkl-2022a

gatk HaplotypeCaller -R $ref \
-I $projectdir${samp}${bamprefix} \
-ploidy 1 \
--emit-ref-confidence BP_RESOLUTION \
-bamout ${samp}.gatk.bam \
-O $samp.g.forbam.vcf

done
else

echo "variants will NOT be called from individual sample BAMs"

fi
#remove the existing list of vcfs

if [ -e "${vcflist}" ]; then
echo "removing the existing list of VCFs $vcflist"
rm $vcflist
fi

#make a list of the gvcfs to compare group genotypes
for samp in $samplist
do
    echo $samp.g.forbam.vcf >> $vcflist
done

#joint genotyping

module purge; module load GATK/4.5.0.0-gimkl-2022a

# first combine gvcfs
gatk CombineGVCFs \
-R $ref \
--variant $vcflist \
-O ${grpname}_all_indivs_variants.vcf

# then do joint genotyping of all indivs within single file
gatk GenotypeGVCFs \
-R $ref \
--variant ${grpname}_all_indivs_variants.vcf \
-O ${grpname}_unfiltered_variants.vcf

#variant recalibration (not recommended when only using mtDNA vars) /filtering

#filter the SNPs
echo "pulling out SNPs and INDELs"

#select snps
gatk SelectVariants \
-R $ref \
--variant ${grpname}_unfiltered_variants.vcf \
--select-type-to-include SNP \
-O raw_snps.vcf 

#select the indels
gatk SelectVariants \
-R $ref \
--variant ${grpname}_unfiltered_variants.vcf \
--select-type-to-include INDEL \
-O raw_indels.vcf 

echo "filtering the SNPs and INDELs"

if [ $filtertype = "hard" ]; then

#filter snps
gatk VariantFiltration \
-R $ref \
--variant raw_snps.vcf \
--filter-expression "QD < 2.0 || FS > 60.0 || MQ < $MQthreshold" \
--filter-name "gatk_snp_filter" \
-O ${grpname}_filtered_snps.vcf 

#filter indels
gatk VariantFiltration \
-R $ref \
-V raw_indels.vcf \
--filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
--filter-name "gatk_indel_filter" \
-O ${grpname}_filtered_indels.vcf 

#recombine indels + snps # GATK's CombineVariants no longer exists, but is available as part of DISCVRSeq
ml purge; ml Java/20.0.2
java -jar /nesi/project/ga03186/bin/DISCVRSeq-1.3.76.jar MergeVcfsAndGenotypes \
-R $ref \
--variant:snp ${grpname}_filtered_snps.vcf \
--variant:ind ${grpname}_filtered_indels.vcf \
-genotypeMergeOption PRIORITIZE \
-priority snp,ind \
-O ${grpname}_samples_include-filtered.vcf

#remove variants that have failed filtering
ml purge; ml GATK/4.5.0.0-gimkl-2022a
gatk SelectVariants \
-R $ref \
--variant ${grpname}_samples_include-filtered.vcf \
--exclude-filtered \
-O ${grpname}_samples.vcf

else
echo "currently not set up to do VQSR filtering, please perform manually"
fi

####################
else
#create bamlist
if [ ! -e "$projectdir$grpname.bam.list" ]; then
echo "Could not find bamlist! Creating one from the samplist"
ls ${projectdir}*${bamprefix} > $projectdir$grpname.TEST.bam.list
else
echo "Found bamlist for group $grpname"
fi

#call the variants 
if [ ! -e "$projectdir${grpname}_samples.vcf" ]; then
echo "No VCF of called variants for this group! Calling a new VCF with HaplotypeCaller. This make take some time"
module purge; module load GATK/4.5.0.0-gimkl-2022a
gatk HaplotypeCaller \
-R $ref \
-I ${projectdir}${bamlist} \
-ploidy $ploidy \
-O $projectdir${grpcalldir}${grpname}_samples.vcf
else
echo "VCF of called variants exists for this group. Please delete the existing VCF if you wish to overwrite"
fi
fi

cd ..
####################
