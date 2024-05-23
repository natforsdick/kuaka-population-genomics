#!/bin/bash -e

# merge the samples sampled in duplicate so replicate fastqs are combined 
cd /nesi/nobackup/ga03186/kuaka-pop-gen/output//03-merged/

cat EXT041-03_S3_val_1.fq.gz D206935_S38_val_1.fq.gz > D206935_S38_EXT041-03_S3_val_1.fq.gz
cat EXT041-03_S3_val_2.fq.gz D206935_S38_val_2.fq.gz > D206935_S38_EXT041-03_S3_val_2.fq.gz
cat EXT041-05_S5_val_1.fq.gz D206949_S60_val_1.fq.gz > D206949_S60_EXT041-05_S5_val_1.fq.gz
cat EXT041-05_S5_val_2.fq.gz D206949_S60_val_2.fq.gz > D206949_S60_EXT041-05_S5_val_2.fq.gz

# and tidy up by archiving the unmerged files so they don't disrupt the pipeline
mkdir non-merged/
mv ./{EXT041-03_S3_val_*.fq.gz,EXT041-05_S5_val_*.fq.gz,D206935_S38_val_*.fq.gz,D206949_S60_val_*.fq.gz} non-merged/
