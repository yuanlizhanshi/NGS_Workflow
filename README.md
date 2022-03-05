# NGS workflow
This is a series of Next generation sequencing (illumina short reads sequencing) workflow created by Snakemake


<!-- TOC titleSize:2 tabSpaces:2 depthFrom:2 depthTo:6 withLinks:1 updateOnSave:1 orderedList:0 skip:0 title:1 charForUnorderedList:* -->
## Table of Contents
* [Setup](#setup)
* [Bulk RNA-seq workflow](#bulk-rna-seq-workflow)
  * [Usage](#usage)
* [Whole genome sequencing workflow](#whole-genome-sequencing-workflow)
<!-- /TOC -->

## Setup
All Requirements software was shown in environment.yaml. Your could install via:
> conda install --file environment.yaml

## Bulk RNA-seq workflow
This is a RNA-seq workflow

It provide two forms of RNA-seq workflow (mapped by STAR or Hisat2).

The diagram of workflow was shown in folder Flow diagram

Before the RNA-seq, the genome index should be built **first**.

The bulid method was shown:

**STAR**:

>STAR \
--runMode genomeGenerate \
--genomeDir ~/STARindex \
--runThreadN 10 \
--genomeFastaFiles chr1.fa \
--sjdbGTFfile chr1.gtf \
--sjdbOverhang 149

**Hisat2**
>hisat2_extract_splice_sites.py chr1.gtf > splice.txt \
>hisat2_extract_exons.py chr1.gtf > exons.txt \
>hisat2-build -p 40  chr1.fa --ss splice.txt --exon exons.txt chr1_index
### Usage
**For dry use**
>snakemake -s RNA-seq.smk -np

**Run workflow**
>snakemake -s RNA-seq.smk -c 4

-c: Default use 4 cores

-----

## Whole genome sequencing workflow

This is a Whole genome sequencing workflow


If you install GATK with conda, please check the version of gatk via:
>gatk --version

if your GATK version is **lower than 4.0**, please install gatk manually from github:
https://github.com/broadinstitute/gatk/releases

note: GTAK rely on java environment, if you don't have java, should use install java first:
>conda instll openjdk -y

Before the WGS, the genome index should be built **first**.

**Bwa**
>bwa index genome.fa

**GATK**
>gatk CreateSequenceDictionary -R genome.fa -O genome.dict
