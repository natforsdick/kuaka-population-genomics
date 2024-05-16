#!/bin/bash -e

#SBATCH -A ga03186
#SBATCH -J conversion
#SBATCH --cpus-per-task=2
#SBATCH --mem 8G
#SBATCH --time=00:20:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=forsdickn@landcareresearch.co.nz
#SBATCH --output %x.%j.out
#SBATCH --error %x.%j.err

cd /nesi/nobackup/ga03186/kuaka-pop-gen/output/05-variant-calling-b/filter-trial/
mkdir -p ./export/

# make bcflist.txt
# ls LD-filter/*.bcf > bcflist.txt
# ls noLD/*.bcf >> bcflist.txt

ml purge
ml BCFtools/1.10.2-GCC-9.2.0 Stacks/2.61-gimkl-2022a

for bcf in $(cat bcflist.txt)
do
    echo bcftools and stacks conversion of $bcf
    base=${bcf%.bcf}
    base=${base##*/}
    bcftools convert -O v -o ./export/${base}.vcf $bcf
    populations -V ./export/${base}.vcf -O ./export/ -M popmap.txt --plink
done

ml purge
ml PLINK/1.09b6.16
for bcf in $(cat bcflist.txt)
do
    echo plink conversion of $bcf
    base=${bcf%.bcf}
    base=${base##*/}
    plink --file ./export/${base}.p.plink --aec --recode A
done
