# Kuaka population genomics

Reference-based population genomics for the critically endangered Kuaka Whenua Hou | Whenua Hou diving petrel (*Pelecanoides whenuahouensis*) and its more common and widespread congener, Kuaka | Common diving petrel (*Pelecanoides urinatrix*). Aims of the project are to assess genetic diversity and differentiation for the two populations, and determine whether observed interspecific hybridisation is having any impact on population recovery of Kuaka Whenua Hou.

This project is a collaboration between the [Department of Conservation](https://www.doc.govt.nz/), [Manaaki Whenua - Landcare Research](https://www.landcareresearch.co.nz/), and [Genomics Aotearoa](https://www.genomics-aotearoa.org.nz/), and reports back to Kāi Tahu – Mana Whenua and Mana Moana (people of the land and people of the sea), Whenua Hou Komiti, and Kaitiaki Rōpū Ki Murihiku.

Collaborators: Nat Forsdick, Johannes Fischer, Igor Debski, Thomas Buckley.

Scripts were originally run on the [NeSI](https://www.nesi.org.nz/) platform via SLURM workload manager, except for R scripts which were run locally. 

This workflow moves through 1) quality control and 2) mapping of individual short-read sequence data, before 3) variant calling against a reference genome for Kuaka Whenua Hou. Variants are filtered, and the resulting SNP sets are exported for 4) analysis and visualisation via genetics packages such as SNPRelate in R. 

## Requirements

* Short-read sequence data for each sample in FASTQ format
* Reference genome assembly in FASTA format

## Software

* [FastQC](https://github.com/s-andrews/FastQC) v0.12.1
* [MultiQC](https://multiqc.info/) v1.13
* [TrimGalore](https://github.com/FelixKrueger/TrimGalore) v0.6.7
* [seqtk](https://github.com/lh3/seqtk) v1.3
* [BWA](https://github.com/lh3/bwa) v0.7.17
* [SAMtools](https://github.com/samtools/samtools) v1.10
* [split_bamfiles_tasks.pl](https://github.com/Lanilen/SubSampler_SNPcaller/split_bamfiles_tasks.pl)
* [BCFtools](https://samtools.github.io/bcftools/bcftools.html) v1.15.1 (or 1.10.2 or 1.19)
* [VCFtools](https://vcftools.sourceforge.net/) v0.1.15
* [Stacks](https://catchenlab.life.illinois.edu/stacks/comp/populations.php) v2.61
* [PLINK](https://www.cog-genomics.org/plink/) v1.09b6.16
* R packages: [SNPRelate](https://github.com/zhengxwen/SNPRelate), vcfR, [adegenet](https://cran.r-project.org/web/packages/adegenet/index.html), [pophelper](https://github.com/royfrancis/pophelper) 
