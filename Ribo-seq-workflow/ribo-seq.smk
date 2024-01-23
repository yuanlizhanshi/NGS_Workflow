from pathlib import Path
import re
import pandas as pd

adapt_seq = 'CTGTAGGCACCATCAAT'


trRNA_index = './mouse_ribo_utils/mm10_trRNA'
RNA_index = './mouse_ribo_utils/mm10_longest_transcript'
ribo_density_script = "./script/calculate_ribosome_density.py"

all_file = Path("rawdata/").glob('*fastq.gz')
all_file_list = [i for i in all_file]

SAMPLES = []
for i in all_file_list:
    file = re.search(r'SRR\d+',i.stem).group()
    if file not in SAMPLES:
        SAMPLES.append(file)
rule all:
    input:
        expand("rawdata/{sample}.fastq.gz",sample=SAMPLES),
        expand("clean_fastq/{sample}_clean.fastq.gz",sample=SAMPLES),
        expand("remove_trRNA/{sample}_trRNA.fastq.gz",sample=SAMPLES),
        expand("ribo_density/{sample}_ribo_density.txt",sample=SAMPLES),
        expand("sortedbam/{sample}.bam",sample=SAMPLES),
        expand("bigwig/{sample}.bigwig",sample=SAMPLES)

rule QC:
    input:
        "rawdata/{sample}.fastq.gz"
    output:
        "clean_fastq/{sample}_clean.fastq.gz",
    threads: 10
    log:
        "clean_fastq/{sample}.html"
    shell:
        "fastp -a {adapt_seq} -w {threads} -i {input} -o {output} --html {log}"

rule remove_trRNA:
    input:
        "clean_fastq/{sample}_clean.fastq.gz"
    output:
        no_trRNA_fq = 'remove_trRNA/{sample}_trRNA.fastq',
        sam = temp('remove_trRNA/{sample}_trRNA.sam')
    threads: 10
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
    threads: 10
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
        "python {ribo_density_script} {input} {output}"


rule samtools_sort:
    input:
        'sam/{sample}.sam'
    output:
        'sortedbam/{sample}.bam'
    threads: 10
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
    threads: 10
    shell:
        'bamCoverage -p {threads} -bs 1 -b {input.bam} -o {output} --normalizeUsing CPM 2>{log}'