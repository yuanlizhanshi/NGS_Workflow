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
* [RNA Chromatin Immunoprecipitation (RIP) sequencing workflow](#rna-chromatin-immunoprecipitation-rip-sequencing-workflow)
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

If you don't have a high-quality reference genome for your species or want to identify more transcripts such as long non-coding RNA, using stringtie after mapping is suggested.



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
This is a Whole genome bisulfite sequencing workflow (from fastq.gz to bed and other downsteam formats).\
There are **two** types of workflow: map by bwa or by bismark, you can choose one of the workflow.


For **bwa-meth**:\
bwa-meth is available at https://github.com/brentp/bwa-meth, please insatall it and its dependencies first\
Then bulid index:
>bwameth.py index genome.fa

The final output file has two formats: methylKit and cytosine_report. The detailed description of the output file format were in https://github.com/dpryan79/MethylDackel

For **bismark**:\
bismark is available at conda and https://github.com/FelixKrueger/Bismark, Genome index also need to build first:
>bismark_genome_preparation ./genome/

Before running this snakemake, there are many parameters that need to change in the snakemake file, such as the aligner, the mismatch number, and whether the strand specific library, etc. you should read the document of bismark carefully.

-----
## Chromatin Immunoprecipitation (ChIP) sequencing workflow

This is a ChIP-seq workflow (from fastq.gz to peak).
Before the ChIP-seq workflow, the genome index of bowtie2 should be built **first**.

>bowtie2-build genome.fa  genome

There are also many parameters need change in snakemake file. Its depends on your ChIP experiment reseach object (Transcription factors or Histone modification), reseach species (Genome size) and replication numbers. you should read the paper of [MACS2](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2008-9-9-r137) carefully.

[Deeptools](https://github.com/deeptools/deepTools) was recommended to generate the meta-plot and heatmap on a reference point or certain genomic regions.

[HOMER](http://homer.salk.edu/homer/ngs/peakMotifs.html) and [MEME suite](http://meme.ebi.edu.au/meme/index.html) was recommended to TF Motif enrichment and analysis.

[ChromHMM](http://compbio.mit.edu/ChromHMM/) was was recommended to peroform chromatin state segmentation.

-----
## ATAC-seq workflow
This is a ATAC-seq workflow (from fastq.gz to peak).

The workflow of the ATAC-seq is similar to the ChIP-seq, so the genome index was the same, the only difference is Peak calling (ATAC **without** input sample).

This workflow is also appropriate for other open chromatin sequencing methods such as Dnase-Seq, MNase-Seq, and FAIRE-Seq, just might need to change the parameters of MACS2.

After Peak calling, To get all peaks information and for downsteam analyse, you should run:

For multi sample, it is recommended to merge peaks with iterative Overlap Peak Merging Procedure, which was first introduced in [Corces & Granja et. al. Science 2018](https://science.sciencemag.org/content/362/6413/eaav1898) and recommonded in [Grandi et al. Nat Protoc 2022](https://www.nature.com/articles/s41596-022-00692-9), the code was modified from [ArchR](https://github.com/GreenleafLab/ArchR).

The differential peaks could calculated with [DEseq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) or [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html). (The identification of differential peaks is still controversial and you should read the literature carefully,such as [DiffBind](http://www.bioconductor.org/packages/release/bioc/html/DiffBind.html)).\
[HINT](http://www.regulatory-genomics.org/hint/introduction/) was recommended to find the TF footprints


-----
## RNA Immunoprecipitation (RIP) sequencing workflow
This is a RIP-seq workflow (from fastq.gz to peak), which is also suit for m6a-seq(meRIP)

The mapping method was consistent with RNA-seq and the Call peak method was consistent with to ChIP-seq. \
For identification of differential peaks, [exomePeak2](https://bioconductor.org/packages/release/bioc/html/exomePeak2.html) was recommended.

-----
## Single cell multi-omics workflow
This is workflow is used for upsteam analysis of Single cell RNA-seq and Single cell ATAC-seq.
Currently only support the **10X genomics** sequencing platform.\
Before the analyse, you should download **Cell Ranger** form the website of [10X genomics](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger).Cell Ranger was an integrated software for Single cell RNA-seq and Single-cell ATAT-seq.
If your used public data, you should **rename** you fastq following the naming rules of cellranger.\
Then index your genomne:
>cellranger mkref \
--genome genome \
--fasta genome.fa \
--genes genome_annotation.gtf \
--nthreads=40

For scRNA-seq:
>cellranger count \
--id sample_ID \
--sample sample_ID \
--localcores 40 \
--no-bam \
--nosecondary \
--fastqs fastq_folder/ \
--transcriptome genome_index

For scATAC-seq:
>cellranger-atac count \
--id=sample_name \
--reference=genome_index \
--fastqs=fastq_folder/ \
--localcores=40

[Seurat](https://satijalab.org/seurat/index.html), [monocell](https://cole-trapnell-lab.github.io/monocle3/), and [ArchR](https://www.archrproject.com/) was recommended to perform downstream analysis of singlecell multi-omics data. 