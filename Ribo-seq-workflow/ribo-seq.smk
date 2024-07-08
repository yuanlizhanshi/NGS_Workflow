from pathlib import Path
import re
import pandas as pd

trRNA_index = '/data/kyh/our_NB4/NB4_ribo/hg38_ribo_utils/hg38_trRNA'
cds_info = './hg38_ribo_utils/cds_info.bed'
UTR5_info = './hg38_ribo_utils/utr5_info.bed'
UTR3_info = './hg38_ribo_utils/utr3_info.bed'
RNA_index = '/data/kyh/our_NB4/NB4_ribo/hg38_ribo_utils/hg38_longest_transcript'
ribo_density_script = "/data/kyh/our_NB4/NB4_ribo/script/calculate_ribosome_density.py"
genome_rna_index = '/data/kyh/our_NB4/NB4_UMI_test/index/hg38_emcc'
longest_transcripts = '/data/kyh/our_NB4/NB4_ribo/hg38_ribo_utils/hg38_longest.transcripts.info.txt'
ultis_path = '/data/kyh/our_NB4/NB4_ribo/hg38_ribo_utils'
gtf = "~/Desktop/hg38/hg38_transcript/genes.gtf"
STAR_index = '~/Desktop/hg38/STAR_ribo_index'
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
        expand("final_bam/{sample}_final.bam",sample=SAMPLES),
        expand("final_bam/{sample}_final.bam.stat",sample=SAMPLES),
        expand("genome_sam/{sample}_STAR_counts.txt",sample=SAMPLES),
        expand("final_bam/{sample}_CDS_counts.bed",sample=SAMPLES),
        expand("final_bam/{sample}_UTR5_counts.bed",sample=SAMPLES),
        expand("final_bam/{sample}_UTR3_counts.bed",sample=SAMPLES)




rule QC:
    input:
        "cutadapt_fastq/Ribo-{sample}.cutadapt.fastq.gz"
    output:
        "clean_fastq/{sample}_clean.fastq.gz",
    threads: 10
    log:
        "clean_fastq/{sample}.html"
    shell:
        "fastp -A -l 20 --length_limit 35 -w {threads} -i {input} -o {output} --html {log}"

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
#here

rule STAR_map:
	input:
		clean_R1 = "remove_trRNA/{sample}_trRNA.fastq.gz"
	output:
		"sam/{sample}/Aligned.sortedByCoord.out.bam"
	threads: 10
	shell:
		"STAR --genomeDir {STAR_index} --twopassMode Basic "
		"--runThreadN {threads} --readFilesCommand zcat "
		"--readFilesIn {input.clean_R1} "
		"--alignEndsType Local "
		"--outSAMtype BAM SortedByCoordinate --sjdbOverhang 30 "
		"--outFilterMultimapNmax 1 --outFilterMismatchNmax 0 "
		"--seedSearchStartLmax 14 --alignIntronMax 10000 "
		"--outFilterIntronMotifs RemoveNoncanonical "
		"--outFileNamePrefix ./sam/{wildcards.sample}/"

rule rename_bam:
	input:
		bam = "sam/{sample}/Aligned.sortedByCoord.out.bam"
	output:
		"sam/{sample}.bam"
	shell:
		"mv {input.bam} sam/{wildcards.sample}.bam"

rule index_STAR_bam:
	input:
		bam = "sam/{sample}.bam"
	output:
		"sam/{sample}.bam.bai"
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

rule count_STAR_bam:
    input:
        bam = 'sam/{sample}.bam',
        bai = 'sam/{sample}.bam.bai'
    output:
        'genome_sam/{sample}_STAR_counts.txt'
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


rule final_bam:
    input:
        'sortedbam/{sample}.bam'
    output:
        'final_bam/{sample}_final.bam'
    shell:
        "samtools view -h -q 10 -b {input} > {output}"

rule inedx_final_bam:
    input:
        'final_bam/{sample}_final.bam'
    output:
        'final_bam/{sample}_final.bam.bai'
    shell:
        "samtools index {input}"

rule stat_final_bam:
    input:
        bam = 'final_bam/{sample}_final.bam',
        bai = 'final_bam/{sample}_final.bam.bai'
    output:
        'final_bam/{sample}_final.bam.stat'
    shell:
        "samtools idxstats {input.bam} > {output}"

rule get_CDS_counts:
    input:
        bam = 'final_bam/{sample}_final.bam',
        bai = 'final_bam/{sample}_final.bam.bai'
    output:
        'final_bam/{sample}_CDS_counts.bed'
    shell:
        "bedtools multicov -bams {input.bam} -bed {cds_info} > {output}"

rule get_utr5_counts:
    input:
        bam = 'final_bam/{sample}_final.bam',
        bai = 'final_bam/{sample}_final.bam.bai'
    output:
        'final_bam/{sample}_UTR5_counts.bed'
    shell:
        "bedtools multicov -bams {input.bam} -bed {UTR5_info} > {output}"

rule get_utr3_counts:
    input:
        bam = 'final_bam/{sample}_final.bam',
        bai = 'final_bam/{sample}_final.bam.bai'
    output:
        'final_bam/{sample}_UTR3_counts.bed'
    shell:
        "bedtools multicov -bams {input.bam} -bed {UTR3_info} > {output}"



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


        
        
