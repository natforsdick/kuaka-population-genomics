#!/bin/bash -e

# prepping the reference file

####################
# variables
# reference file
refdir=/nesi/nobackup/ga03186/kuaka-genome/mitohifi/
reffile="kuaka_final_mitogenome-80bprem.fasta NC_052809.1-P-urinatrix-mitogenome.fasta NC_052809.1-P-urinatrix-mitogenome-50bprem.fasta"

# software
module purge
module load picard/2.26.10-Java-11.0.4 BWA/0.7.18-GCC-12.3.0 

cd $refdir

for ref in $reffile
do
####################
#index the reference fasta file
if [ ! -e $ref.amb ]; then
echo "Index file of reference does not exist: creating index with BWA"
bwa index $ref
else
echo "BWA Index file found"
fi

#Create a Sequence Dictionary if necessary
if [ ! -e $ref.dict ]; then
echo "SequenceDictionary of reference does not exist: creating one with picard"
picard CreateSequenceDictionary \
-R $ref \
-O $ref.dict
else
echo "SequenceDictionary found"
fi
done
