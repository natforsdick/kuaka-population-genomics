# Mapping 

In this section, we map the trimmed, merged fastq data against the reference genome. We then collect mapping statistics. Both mapping and collecting stats are run via arrays.

* `reffile`: reference genome in FASTA format
* `refdir`: location of reference genome
* `fq_list`: list of input trimmed merged data for each sample found in the input directory `INDIR` in FASTQ format


## Example `fq_list` file

```
./sampleID1_val_1.fq.gz
./sampleID2_val_1.fq.gz
./sampleID3_val_1.fq.gz
...
```

## Example `BAMLIST` for `03-merge-stats.sl`

```
/path/to/04-mapped/bam/sampleID1.aligned.sorted.bam
/path/to/04-mapped/bam/sampleID2.aligned.sorted.bam
/path/to/04-mapped/bam/sampleID3.aligned.sorted.bam
...
```
