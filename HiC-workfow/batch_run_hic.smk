import glob
import os
import re

directory_path = "fastq"
samples = []
files = glob.glob(directory_path + '/*.fastq.gz')


for i in files:
    filename = re.search('[A-Z0-9]+',os.path.splitext(os.path.splitext(os.path.basename(i))[0])[0]).group()
    samples.append(filename)

genome_index = '/home/kyh/Desktop/hg38/hg38_bwa/hg38'

digest_bed = '/data/kyh/islet_regulatome/HiC/hg38_digest2.bed'

rule all:
    input:
        expand("fastq/{sample}_1.fastq.gz",sample=samples),
        expand("fastq/{sample}_2.fastq.gz",sample=samples),
        expand("valid_pairs/{sample}.valid_pairs.txt.gz",sample=samples),
        expand("hic/{sample}.hic",sample=samples)
        
        
rule QC:
    input:
        raw_R1 = "fastq/{sample}_1.fastq.gz",
        raw_R2 = "fastq/{sample}_2.fastq.gz"
    output:
        clean_R1 = "clean_fastq/{sample}_1.fastq.gz",
        clean_R2 = "clean_fastq/{sample}_2.fastq.gz"
    log:
        html = "clean_fastq/{sample}.html",
        json = "clean_fastq/{sample}.json",
    threads: 10
    shell:
        "fastp -w {threads} -i {input.raw_R1} -o {output.clean_R1} "
        "-I {input.raw_R2} -O {output.clean_R2} --detect_adapter_for_pe --html {log.html} --json {log.json}"
        
rule mapping:
    input:
        R1 = 'clean_fastq/{sample}_1.fastq.gz',
        R2 = 'clean_fastq/{sample}_2.fastq.gz'
    output:
        'valid_pairs/{sample}.valid_pairs.txt'
    threads: 20
    shell:
        " ./script/bwa_map2.sh -p {threads} -1 {input.R1} -2 {input.R2} -g {genome_index} -e {digest_bed} -n {wildcards.sample} "

rule gzip:
    input:
        'valid_pairs/{sample}.valid_pairs.txt'
    output:
        'valid_pairs/{sample}.valid_pairs.txt.gz'
    threads: 20
    shell:
        "bgzip -@ {threads} {input}"
rule filter_by_mapq:
    input:
        'valid_pairs/{sample}.valid_pairs.txt.gz'
    output:
        'valid_pairs/{sample}.valid_pairs_mapq30.txt.gz'  
    threads: 20
    shell:
        "zcat {input} |awk ' $10 >= 30 && $11 >= 30 {{print $0}}' |bgzip -@ {threads} > {output}"

rule pairs_dedup:
    input:
        'valid_pairs/{sample}.valid_pairs_mapq30.txt.gz'
    output:
        'valid_pairs/{sample}.valid_pairs_mapq30_dedup_sort.txt.gz'  
    threads: 20
    shell:
        "zcat {input} | ./script/sort -T ./ -k3,3 -k7,7 -k4,4n -k8,8n -u | bgzip -@ {threads} > {output} "
       
rule pairs2hic:        
    input:
        'valid_pairs/{sample}.valid_pairs_mapq30_dedup_sort.txt.gz'
    output:
        'hic/{sample}.hic'
    threads: 20
    shell:
        "java -Xmx1000G -jar ./script/juicer_tools_1.22.01.jar pre --threads {threads} -q 30 -f {digest_bed} {input} {output} hg38"
