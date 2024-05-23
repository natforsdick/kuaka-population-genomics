# Raw data QC and pre-processing

These scripts are used to QC the raw population data, perform quality trimming, and merge data where samples were sequenced on multiple Illumina lanes. A small number of samples that produced disproportionately large amounts of data are downsampled, and a data for a small number of samples that were sequenced in duplicate are merged.

The spread of sequence depth led us to partition the sample data into 'large', 'mid', and 'small' data sets, to facilitate more efficient use of SLURM resources.

Example format of sample lists used as input here:

```
# 01-fastqc-large-array.sl and 02-trimgalore.sl - all samples, split into small/mid/large-fastq.txt
sequenceID1-sampleID1_L001_R1_001.fastq.gz
sequenceID1-sampleID1_L002_R1_001.fastq.gz
sequenceID2-sampleID2_L001_R1_001.fastq.gz
sequenceID2-sampleID2_L002_R1_001.fastq.gz
...

# 04-subsamp.sl - the outputs of trimgalore for those samples that produced substantially higher sequencing depth than the mean, in subsamp-in.txt
sampleID1_val_1.fq.gz
sampleID2_val_1.fq.gz
sampleID3_val_1.fq.gz
...
```
