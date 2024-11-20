#!/bin/bash -e

# call FASTA files from BAM files produced by 03-variant-calling-mito.sl
# Uses GATK's HaploCaller to call best variants at a set depth
# Nullifies any bases (ref or var) at or less than depth specified by the 'mincoverage' variable

####################
# variables to change 
# location of the FASTA reference file
refdir=/nesi/nobackup/ga03186/kuaka-genome/mitohifi/
reffile=kuaka_final_mitogenome-80bprem.fasta
ref=$refdir$reffile
#name of this group of samples
grpname=KWH-mito
#project directory where input vcfs are
projectdir=/nesi/nobackup/ga03186/kuaka-pop-gen/output/mitogenomes/KWH/
#the ploidy of the data (eg, 1 for mitochondria, 2 for autosomal)
ploidy=1
#the MINIMUM coverage allowed to call a base for these samples (eg, if set to 1, all bases covered by 1 or more reads will be in FASTA)
mincoverage="10"
# mitogenome ref sequence name in file
mtchr="D212900-mitogenome"

#below: a list of sample names (in any order) as they were processed by the fastq_pipeline_loop[version].sh
#these sample names currently should be in double quotes (around all) and separated by whitespace (spaces or tabs)
# samplist containing list of file prefixes for the trimmed short-read data
indir=/nesi/nobackup/ga03186/kuaka-pop-gen/output/03-merged/
samplist=${indir}mitolist-KWH-nolow.txt

#variables that shouldn't need changing
bamlist=$grpname.bam.list
#oneunder=$((mincoverage - 1))
gencalldir=/nesi/nobackup/ga03186/kuaka-pop-gen/output/mitogenomes/KWH/call_genotypes-KWH-only/
grpcalldir="call_group_genotypes/"
#specify the name of the BAM file to call from
bamprefix="_remdup.bam"

####################
#enter the project directory

if [[ ! -f /nesi/project/ga03186/kuaka-pop-gen/05-mitogenomes/filter_cover_for_vcfs.R ]] ; then
    echo "Could not find filter_cover_for_vcfs.R in main processing directory, aborting."
    exit
fi

cd $projectdir

#create dictionary file of the reference if it doesn't exist
dictfile=${reffile%.*}.dict
if [ ! -e "$refdir${reffile%.*}.dict" ]; then
    echo "Dictionary file of reference does not exist: creating dictionary"
    module purge; module load picard/2.26.10-Java-11.0.4
    picard CreateSequenceDictionary -R $ref -O $refdir${reffile%.*}.dict
else
    echo "Dictionary file found"
fi

#index the reference if it is not indexed already
idxfile=${reffile%.*}.fai
if [ ! -e "$refdir$reffile.fai" ]; then
    echo "Reference is not indexed: indexing with SAMTOOLS"
    module purge; module load SAMtools/1.19-GCC-12.3.0
    samtools faidx $ref
else
    echo "Index for reference file found"
fi

#call variants across entire group
if [ ! -e "${projectdir}${grpcalldir}${grpname}_samples.vcf" ]; then
    echo "No VCF of called variants found for this group! Please run variantcall_group_vcf_from_bam script first"
else
    echo "Using the called variants from ${grpname}_samples.vcf"
fi

####################
#begin processing each sample
mkdir -p $gencalldir
cd $gencalldir
for samp in $(cat $samplist)
do
    #enter directory
    echo "starting FASTA creation for sample $samp"
    mkdir -p ${gencalldir}${samp}_process
    cd ${gencalldir}${samp}_process

    if [ -e "vcf_list.list" ]; then
    echo "Removing an existing list of VCF coverage files"
    rm ${gencalldir}vcf_list.list
    fi

    #get the VCF for the BEST variant calls
    echo "creating a VCF of $samp from the group VCF"
    module purge; module load GATK/4.5.0.0-gimkl-2022a
    gatk SelectVariants \
        -R $ref \
        --variant ${projectdir}${grpcalldir}${grpname}_samples.vcf \
        --sample-name $samp \
        --exclude-non-variants \
        -O ${samp}_haploAll.vcf


    #make VCF
    echo "looking for the BP_RESOLUTION VCF for $samp"

    if [ -e "${projectdir}${grpcalldir}${samp}.g.forbam.vcf" ]; then
        echo "found BP_RESOLUTION VCF for $samp"
        cp ${projectdir}${grpcalldir}${samp}.g.forbam.vcf .
        BPvcf=${samp}.g.forbam.vcf
    else 
        echo "Cannot find BP_RESOLUTION VCF for $samp! Creating a VCF of the entire reference for $samp"

        gatk HaplotypeCaller \
            -R $ref \
            -I ../${samp}${bamprefix} \
            -ploidy 1 \
            --emit-ref-confidence BP_RESOLUTION \
            -O $samp.v.vcf

        BPvcf=$samp.v.vcf
    fi
    #

    #make a hack vcf from the ERC vcf
    echo "selecting low coverage variants from the VCF for $samp"

    #get the VCF
    vcfheadl=$(awk '/#CHROM/{ print NR; exit }' $BPvcf)

    #get coverage
    gatk DepthOfCoverage \
        -R $ref \
        -O ${samp}_coverage_gatk_table \
        --output-format TABLE \
        -L $mtchr \
        -I ${projectdir}call_group_genotypes/${samp}.gatk.bam

    #get the vcf without the header
    tail -n +$vcfheadl $BPvcf > $samp.nohead.vcf

    #run the VCFs
    for cov in $mincoverage
    do
        #run the script to create a min-cov vcf
        echo "running the min.cov script for coverage $cov"
        module purge; module load R/4.3.1-gimkl-2022a
        Rscript /nesi/project/ga03186/kuaka-pop-gen/05-mitogenomes/filter_cover_for_vcfs.R ${projectdir}../coverage/${samp}_coverage.txt $samp.nohead.vcf $cov

        #recreate the vcf
        BPvcf=${samp}.g.forbam.vcf
        head -n $vcfheadl $BPvcf > $samp.new$cov.vcf
        cat $samp.new$cov.vcf gens.$cov.vcf > $samp.new$cov.full.vcf

        #remove the NON_REFs from the vcfs
        vcffile=$samp.new$cov.full.vcf
        echo "altering $vcffile"
        sed -i 's/,<NON_REF>//g' $vcffile
        sed -i 's/A\t<NON_REF>/N\tA/g' $vcffile
        sed -i 's/C\t<NON_REF>/N\tC/g' $vcffile
        sed -i 's/G\t<NON_REF>/N\tG/g' $vcffile
        sed -i 's/T\t<NON_REF>/N\tT/g' $vcffile

        #make fasta for each coverage
        echo "creating a FASTA of bases above $cov reads for $samp"
        module purge; module load GATK/4.5.0.0-gimkl-2022a
        gatk IndexFeatureFile \
            -I $vcffile

        gatk FastaAlternateReferenceMaker \
            -R $ref \
            --variant ${samp}_haploAll.vcf \
            --snp-mask $vcffile \
            --snp-mask-priority \
            -L $mtchr \
            -O ${samp}.min${cov}.filt.fasta
        done
    cd ${gencalldir}
done
