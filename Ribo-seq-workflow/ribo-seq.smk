from pathlib import Path
import re
import pandas as pd

trRNA_index = 'hg38_trRNA'

RNA_index = 'hg38_longest_transcript'
ribo_density_script = "script/calculate_ribosome_density.py"
genome_rna_index = 'hg38_emcc'
longest_transcripts = 'hg38_longest.transcripts.info.txt'
ultis_path = 'hg38_ribo_utils'
gtf = "~/Desktop/hg38/hg38_transcript/genes.gtf"

all_file = Path("cutadapt_fastq/").glob('*fastq.gz')
all_file_list = [i for i in all_file]

SAMPLES = []
for i in all_file_list:
    file = re.search(r'(?<=-)\w+_\d',i.stem).group()
    if file not in SAMPLES:
        SAMPLES.append(file)
        
        
rule all:
    input:
        expand("cutadapt_fastq/Ribo-{sample}.cutadapt.fastq.gz",sample=SAMPLES),
        expand("clean_fastq/{sample}_clean.fastq.gz",sample=SAMPLES),
        expand("remove_trRNA/{sample}_trRNA.fastq.gz",sample=SAMPLES),
        expand("ribo_density/{sample}_ribo_density.txt",sample=SAMPLES),
        expand("sortedbam/{sample}.bam",sample=SAMPLES),
        expand("genome_sam/{sample}.bam",sample=SAMPLES),
        expand("genome_sam/{sample}.bam.bai",sample=SAMPLES),
        expand("reads_length/{sample}_reads_length.txt.gz",sample=SAMPLES),
        expand("reads_length/{sample}_start_periodicity.pdf",sample=SAMPLES),
        expand("genome_sam/{sample}_counts.txt",sample=SAMPLES),
        expand("bigwig/{sample}.bigwig",sample=SAMPLES)

rule QC:
    input:
        "cutadapt_fastq/Ribo-{sample}.cutadapt.fastq.gz"
    output:
        "clean_fastq/{sample}_clean.fastq.gz",
    threads: 10
    log:
        "clean_fastq/{sample}.html"
    shell:
        "fastp -A -l 25 --length_limit 35 -w {threads} -i {input} -o {output} --html {log}"

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

rule gzip_reads:
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


rule genome_map:
    input:
        "remove_trRNA/{sample}_trRNA.fastq.gz"
    output:
        sam = temp('genome_sam/{sample}.sam')
    log:
        "genome_sam/{sample}_genome_mapping.log"
    threads: 10
    shell:
        "hisat2 --threads {threads} -x {genome_rna_index} -U {input} -S {output.sam} 2> {log}"
        
rule sort_genome_bam:
    input:
        'genome_sam/{sample}.sam'
    output:
        'genome_sam/{sample}.bam'
    threads: 20
    shell:
        'samtools sort -@ {threads} -o {output} {input}'


rule index_genome_bam:
    input:
        'genome_sam/{sample}.bam'
    output:
        'genome_sam/{sample}.bam.bai'
    shell:
        "samtools index {input}"
 
rule count_genome_bam:
    input:
        bam = 'genome_sam/{sample}.bam',
        bai = 'genome_sam/{sample}.bam.bai'
    output:
        'genome_sam/{sample}_counts.txt'
    threads: 1
    shell:
        'ModifyHTseq -i {input.bam} -g {gtf} -o {output} -t exon -m union -q 10 --minLen 25 --maxLen 35 --exclude-first 45 --exclude-last 15 --id-type gene_name'


        
    
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

rule plot_reads:
    input:
        bam = 'sortedbam/{sample}.bam',
        bai =  'sortedbam/{sample}.bam.bai'
    output:
        "reads_length/{sample}_reads_length.txt"
    params:
        output_prefix = "reads_length/{sample}"
    shell:
        "LengthDistribution -i {input.bam} -o {params.output_prefix} -f bam"

rule gzip:
    input:
        "reads_length/{sample}_reads_length.txt"
    output:
        "reads_length/{sample}_reads_length.txt.gz"
    shell:
        "gzip {input} {output}"


rule plot_Periodicity:
    input:
        bam = 'sortedbam/{sample}.bam',
        bai =  'sortedbam/{sample}.bam.bai'
    output:
        "reads_length/{sample}_start_periodicity.pdf"
    params:
        output_prefix = "reads_length/{sample}"
    shell:
        "Periodicity -i {input.bam} -a {ultis_path} -o {params.output_prefix}  -c {longest_transcripts} -L 28 -R 30"

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
        
        
