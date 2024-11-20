#!/bin/bash -e
#SBATCH -A ga03186
#SBATCH -J admix
#SBATCH --time 07:00:00 #
#SBATCH -c 16
#SBATCH -n 1
#SBATCH --mem=8GB
#SBATCH --array=1-10
#SBATCH --mail-user=forsdickn@landcareresearch.co.nz
#SBATCH --mail-type=FAIL,END
#SBATCH --output admix_4.%j.out # CHANGE number for new run
#SBATCH --error admix_4.%j.err #  CHANGE number for new run

# To run admixture for investigating species delimitation, population structuring, and introgression.
# Make directories for the outputs
mkdir -p /nesi/nobackup/ga03186/kuaka-pop-gen/output/admixture/admixture_1007_${SLURM_ARRAY_TASK_ID}
cd /nesi/nobackup/ga03186/kuaka-pop-gen/output/admixture/admixture_1007_${SLURM_ARRAY_TASK_ID}

# Call admixture program
admixture=/nesi/project/ga03186/bin/admixture_linux-1.3.0/admixture ## to update
INDIR=/nesi/nobackup/ga03186/kuaka-pop-gen/output/05-variant-calling-b/filter-trial-b/export/
INBED=Petrel_VariantCalls_5x_coverage_0.1site_missing_noMinGQ.bcf.recode_0.6LD_VariantCalls.bed
filename=$(basename $INBED .bed)
# Designate how many hypothetical populations/groups to test:
for k in {1..10};
do

# Run admixture. -C = designate termination criteria, --cv = 10-fold cross-validation. 
# -s = seed, -j = number of threads (multithreaded mode).
$admixture -C 0.0001 --cv=10 ${INDIR}${INBED} $k \
-s time -j24 | tee log_${filename}_${k}.out 

done
