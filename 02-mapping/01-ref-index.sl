#!/bin/bash -e
#SBATCH -J index
#SBATCH -A ga03186
#SBATCH --time=00:30:00 
#SBATCH --mem=4G 
#SBATCH --cpus-per-task=2
#SBATCH --out=%x.%j.out
#SBATCH --err=%x.%j.err

reffile=01-kuaka-hifiasm-p_ctg-purged-clean-omnic-mapped.PT-yahsNMC_scaffolds_final
refdir=/nesi/nobackup/ga03186/kuaka-genome/05-scaffolding/05b-Dovetail-OmniC/all-data-yahs/
REF=$refdir$reffile

###########

module purge
module load  BWA/0.7.17-GCC-9.2.0

###########

# Index reference
if [ ! -f ${REF}.amb  ]; then
    cd $refdir
    echo "The reference has not been indexed. Indexing now"

    bwa index -a bwtsw $REF.fa
    else
        echo "BWA index file found" 
fi

echo correcting output suffixes
mv ${REF}.fa.sa ${REF}.sa
mv ${REF}.fa.amb ${REF}.amb
mv ${REF}.fa.ann ${REF}.ann
mv ${REF}.fa.pac ${REF}.pac
mv ${REF}.fa.bwt ${REF}.bwt
ls -lt $refdir
