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

## Notes on filtering

I recommend that prior to running any filtering, that you get a sense for the variants returned so you can set reasonable thresholds. The two key parameters to look at here ae raw read deph (DP) and average mapping quality (MQ). Here I am only interested in those sites that are not indels, so we exclude those. To pull these out from our final sorted variant file (in this case `Petrel_VariantCalls_concat_sort.bcf`):

```bash
module load BCFtools/1.19-GCC-11.3.0 # load BCFtools

# get dp
bcftools view -H Petrel_VariantCalls_concat_sort.bcf | grep -v "INDEL" | awk '{print $8}' | cut -d ";" -f1 | cut -d "=" -f2 > DP.txt

# get mq
bcftools view -H Petrel_VariantCalls_concat_sort.bcf | grep -v "INDEL" | sed -n -e 's/^.*MQ=//p' | cut -f1 > MQ.txt

# combine into a single output file
paste DP.txt MQ.txt > Petrel_VariantCalls_concat_sort-DP-MQ.txt

# then add DP and MQ as header using your preferred text editor
```

Then import into R and generate some stats:

```r
rawstats <- read.delim("./Petrel_VariantCalls_concat_sort-DP-MQ.txt")
head(rawstats)

summary(rawstats)

# plot DP distribution

w <- rawstats$DP

h<-hist(w, breaks=20, col="red", xlab="Raw read depth per variant",
   main="DP: Histogram with Normal Curve")
xfit<-seq(min(w),max(w),length=40)
yfit<-dnorm(xfit,mean=mean(w),sd=sd(w))
yfit <- yfit*diff(h$mids[1:2])*length(w)
lines(xfit, yfit, col="blue", lwd=2)

# plot MQ distribution
x <- rawstats$MQ

h<-hist(x, breaks=10, col="red", xlab="Average mapping quality per variant",
   main="MQ: Histogram with Normal Curve")
xfit<-seq(min(x),max(x),length=40)
yfit<-dnorm(xfit,mean=mean(x),sd=sd(x))
yfit <- yfit*diff(h$mids[1:2])*length(x)
lines(xfit, yfit, col="blue", lwd=2)
```

In this case, both are skewed to the high end of the distribution, so we want to do a bit more digging to understand what cutoffs will be appropriate. To do this, I want to know what proportion of sites are lower than the respective means.

```r
# get total number of variants
a <- nrow(rawstats)

# replace 617 with your mean DP value
b <- nrow(rawstats[rawstats$DP<617, ])
# calculate this as a percentage of variants retained
b/a*100

# replace 56 with your mean MQ value
e <- nrow(rawstats[rawstats$MQ<56, ])
# calculate this as a percentage of variants retained
e/a*100
```

You can calculate these for the lower ends of the distribution to determine appropriate filtering cutoffs. For DP, we also want to consider setting a maximum cutoff, to take into account those regions that may have high depth due to collapsed alignments of repetitive regions etc. There are other ways to filter the data too, but it's helpful to have a couple of things that you can lock in to the filtering step ([03-filtering.sh](./03-filtering.sh)). 



