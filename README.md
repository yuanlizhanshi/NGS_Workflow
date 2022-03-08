# NGS workflow
This is a series of Next generation sequencing (illumina short reads sequencing) workflow created by Snakemake


<!-- TOC titleSize:2 tabSpaces:2 depthFrom:2 depthTo:6 withLinks:1 updateOnSave:1 orderedList:0 skip:0 title:1 charForUnorderedList:* -->
## Table of Contents
* [Setup](#setup)
* [Usage](#usage)
* [Bulk RNA Sequencing workflow](#bulk-rna-sequencing-workflow)
* [Whole Genome Sequencing workflow](#whole-genome-sequencing-workflow)
* [Whole Genome Bisulfite Sequencing workflow](#whole-genome-bisulfite-sequencing-workflow)
* [Chromatin Immunoprecipitation (ChIP) sequencing workflow](#chromatin-immunoprecipitation-chip-sequencing-workflow)
* [ATAC-seq workflow](#atac-seq-workflow)
* [Single cell multi-omics workflow](#single-cell-multi-omics-workflow)
<!-- /TOC -->

## Setup
All workflow all have its requirements softwares which were shown in environment.yaml. You could install via:
> conda install --file environment.yaml

## Usage

**For dry use**
>snakemake -s snakemake.smk -np

**Run workflow**
>snakemake -s snakemake.smk -c 4

-c: Default use 4 cores

## Bulk RNA Sequencing workflow
This is a RNA-seq workflow (from fastq.gz to count)

It provide two forms of RNA-seq workflow (mapped by STAR or Hisat2).

The diagram of workflow was shown in folder Flow diagram

Before the RNA-seq, the genome index should be built **first**.

The bulid method was shown:

**STAR**:

>STAR \
--runMode genomeGenerate \
--genomeDir ~/STARindex \
--runThreadN 40 \
--genomeFastaFiles chr1.fa \
--sjdbGTFfile chr1.gtf \
--sjdbOverhang 149

**Hisat2**
>hisat2_extract_splice_sites.py chr1.gtf > splice.txt \
hisat2_extract_exons.py chr1.gtf > exons.txt \
hisat2-build -p 40  chr1.fa --ss splice.txt --exon exons.txt chr1_index

-----
## Whole Genome Sequencing workflow

This is a Whole genome sequencing workflow (from fastq.gz to vcf.gz)\
**Note**: This vcf.gz needs further quality control in at further study.

If you install GATK with conda, please check the version of gatk via:
>gatk --version

if your GATK version is **lower than 4.0**, please install gatk manually from github:https://github.com/broadinstitute/gatk/releases

**Note**: GTAK rely on java environment, if you don't have java, should use install java first:
>conda instll openjdk -y

Before the WGS, the genome index should be built **first** by bwa, samtools,and GATK:
>bwa index genome.fa \
samtools faidx genome.fa \
gatk CreateSequenceDictionary -R genome.fa -O genome.dict

-----
## Whole Genome Bisulfite Sequencing workflow
This is a Whole genome bisulfite sequencing workflow (from fastq.gz to bed and other downsteam formats)\
There are **two** types of workflow: map by bwa or by bismark, you can choose one of the workflow.


For **bwa-meth**:\
bwa-meth is available at https://github.com/brentp/bwa-meth, please insatall it and its dependencies first\
Then bulid index:
>bwameth.py index genome.fa

The final output file has two formats: methylKit and cytosine_report. The detailed description of the output file format were in https://github.com/dpryan79/MethylDackel

For **bismark**:\
bismark is available at conda and https://github.com/FelixKrueger/Bismark, Genome index also need to build first:
>bismark_genome_preparation ./genome/

Before running this snakemake, there are many parameters that need to change in the snakemake file, such as the aligner, the mismatch number, and whether the strand specific library, you should read the document of bismark carefully.

-----
## Chromatin Immunoprecipitation (ChIP) sequencing workflow


-----
## ATAC-seq workflow

-----
## Single cell multi-omics workflow
This is workflow is used for upsteam analysis of Single cell RNA-seq and Single cell ATAC-seq.
Currently only support the **10X genomics** sequencing platform.\
Before the analyse, you should download **Cell Ranger** form the website of [10X genomics](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger).
