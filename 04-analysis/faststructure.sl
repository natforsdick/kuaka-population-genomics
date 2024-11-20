#!/bin/bash -e
#SBATCH -A ga03186
#SBATCH -J faststr_2
#SBATCH --time 3:00:00 #
#SBATCH -c 2
#SBATCH --mem=3G
#SBATCH --mail-user=forsdickn@landcareresearch.co.nz
#SBATCH --mail-type=FAIL,END
#SBATCH --output faststr_1007.%j.out # CHANGE number for new run
#SBATCH --error faststr_1007.%j.err #  CHANGE number for new run
#SBATCH --array=1-9

module load fastStructure/1.0-gimkl-2020a-Python-2.7.18

INDIR=/nesi/nobackup/ga03186/kuaka-pop-gen/output/05-variant-calling-b/filter-trial-b/export/
INBED=Petrel_VariantCalls_5x_coverage_0.1site_missing_noMinGQ.bcf.recode_0.6LD_VariantCalls.bed
filename=$(basename $INBED .bed)

echo "Beginning faststructure run for biallelic set for $k at "$(date)
mkdir -p /nesi/nobackup/ga03186/kuaka-pop-gen/output/faststructure/faststr_1007_${SLURM_ARRAY_TASK_ID}
cd /nesi/nobackup/ga03186/kuaka-pop-gen/output/faststructure/faststr_1007_${SLURM_ARRAY_TASK_ID}

for K in {1..10};
do
    structure.py -K $K --input=${INDIR}${filename} \
        --output=${filename}_${K} \
        --cv=100 --full --prior=simple
done
