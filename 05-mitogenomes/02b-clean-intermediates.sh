#!/bin/bash -e
# cleaning up intermediate files from mitomapping

# path to trimmed short-read data
indir=/nesi/nobackup/ga03186/kuaka-pop-gen/output/03-merged/
# path to mitomapping output directory
outdir=/nesi/nobackup/ga03186/kuaka-pop-gen/output/mitogenomes/
#outdir=/nesi/nobackup/ga03186/kuaka-pop-gen/output/mitogenomes/CDP/

# samplist containing list of file prefixes for the trimmed short-read data
samplist=mitolist.txt
fq1=_1.fq
fq2=_2.fq

cd $indir

for samp in $(cat $samplist)
do
    if [ ! -e $indir${samp}${fq1}.gz ] || [ ! -e $indir${samp}${fq2}.gz ]; then
        echo "Zipping fastq files"
        gzip $indir${samp}${fq1}
        gzip $indir${samp}${fq2}
    fi

    rm ${outdir}${samp}_sort.bam
    rm ${outdir}${samp}_sort.bam.bai
    rm ${outdir}${samp}_allsortclean.bam
    rm ${outdir}${samp}_maponly.bam 
    rm ${outdir}${samp}_allsortclean.bai
done
