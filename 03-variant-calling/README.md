# Variant calling

Here we split mapped outputs into chunks for efficient processing through the BCFtools MPILEUP pipeline to call variants for all samples against the reference genome. We then process the outputs into a single all-sample coordinate-sorted variant file. Variants are then filtered to extract high-quality bi-allelic SNPs. The final SNP sets are converted to various output formats for downstream analyses, and a range of statistics are collected.

## Software

* BCFtools/1.15.1 or 1.10.2 or 1.19
* SAMtools/1.15.1
* VCFtools/0.1.15
* Stacks/2.61
* PLINK/1.09b6.16
* [split_bamfiles_tasks.pl](https://github.com/Lanilen/SubSampler_SNPcaller/split_bamfiles_tasks.pl)

## Inputs

* `ref`: `/path/to/reference-genome.fa`
* `INBCF`: The all-sample BCF output from `02-pre-varfiltering.sl`

## Example list files

### `samplist.txt` for [01-mpileup.sl](01-mpileup.sl) and [02-pre-varfiltering.sl](02-pre-varfiltering.sl)

This file lists the aligned sorted BAM files for each sample produced by [02-readgroup-mapping.sl](../02-mapping/02-readgroup-mapping.sl).

```
/path/to/04-mapped/bam/sampleID1.aligned.sorted.bam
/path/to/04-mapped/bam/sampleID2.aligned.sorted.bam
/path/to/04-mapped/bam/sampleID3.aligned.sorted.bam
...
```

### `bcflist.txt` for [02-pre-varfiltering.sl](02-pre-varfiltering.sl)

This file lists the BCF files for each chunk produced in [01-mpileup.sl](01-mpileup.sl), and is automatically produced during [02-pre-varfiltering.sl](02-pre-varfiltering.sl).  

```
/path/to/05-variant-calling/species_1_raw_VariantCalls.bcf.gz
/path/to/05-variant-calling/species_2_raw_VariantCalls.bcf.gz
/path/to/05-variant-calling/species_3_raw_VariantCalls.bcf.gz
...
```

### `bcflist.txt` for [04-export-format.sl](04-export-format.sl)

This file lists the filtered BCF outputs. To make this file:

```
cd /path/to/output/05-variant-calling/filter-trial/
ls LD-filter/*.bcf > bcflist.txt
ls noLD/*.bcf >> bcflist.txt
```
