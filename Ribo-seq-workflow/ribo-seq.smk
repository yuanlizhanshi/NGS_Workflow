import glob
import os
import re

directory_path = "clean_fastq"
trRNA_index = '/mnt/g/ribo_test/index/mm_trRNA_bt2_index'
RNA_index = '/mnt/g/ribo_test/index/mm_longest_transcript'
fastq = []
files = glob.glob(directory_path + '/*.gz')
for i in files:
    filename = os.path.splitext(os.path.basename(i))[0]
    filename = re.search('SRR\d+',filename).group()
    fastq.append(filename)

rule all:
    input:
        expand("clean_fastq/{sample}_clean.fastq.gz",sample=fastq),
        expand("remove_trRNA/{sample}_trRNA.fastq.gz",sample=fastq),
        expand("ribo_density/{sample}_ribo_density.txt",sample=fastq),
        expand("sortedbam/{sample}.bam",sample=fastq),
        expand("bigwig/{sample}.bigwig",sample=fastq)

rule remove_trRNA:
    input:
        "clean_fastq/{sample}_clean.fastq.gz"
    output:
        no_trRNA_fq = 'remove_trRNA/{sample}_trRNA.fastq',
        sam = temp('remove_trRNA/{sample}_trRNA.sam')
    threads: 8
    log:
        'remove_trRNA/{sample}_trRNA.log'
    shell:
        "bowtie2 --threads {threads} -x {trRNA_index} -U {input} -S {output.sam} --un {output.no_trRNA_fq} 2> {log}"

rule gzip:
    input:
        "remove_trRNA/{sample}_trRNA.fastq"
    output:
        "remove_trRNA/{sample}_trRNA.fastq.gz"
    threads: 1
    shell:
        "gzip {input}"

rule align_RNA:
    input:
        "remove_trRNA/{sample}_trRNA.fastq.gz"
    output:
        sam = temp('sam/{sample}.sam')
    threads: 8
    log:
        'sam/{sample}.log'
    shell:
        "hisat2 --threads {threads} -x {RNA_index} -U {input} -S {output.sam} 2> {log}"


rule calculate_ribo_density:
    input:
        "sam/{sample}.sam"
    output:
        'ribo_density/{sample}_ribo_density.txt'
    threads: 1
    shell:
        "python ./script/calculate_ribosome_density.py {input} {output}"


rule samtools_sort:
    input:
        'sam/{sample}.sam'
    output:
        'sortedbam/{sample}.bam'
    threads: 8
    shell:
        'samtools sort -@ {threads} -o {output} {input}'

rule bam_index:
    input:
        'sortedbam/{sample}.bam'
    output:
        'sortedbam/{sample}.bam.bai'
    shell:
        "samtools index {input}"

rule bam_to_bigwig:
    input:
        bam = 'sortedbam/{sample}.bam',
        bam_index = 'sortedbam/{sample}.bam.bai'
    output:
        'bigwig/{sample}.bigwig'
    log:
        'bigwig/{sample}.bigwig.log'
    shell:
        'bamCoverage -bs 1 -b {input.bam} -o {output} --normalizeUsing CPM 2>{log}'
