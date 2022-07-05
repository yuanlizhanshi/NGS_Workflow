SAMPLES = {'test_1','test_2'}
genome = '/home/kyh/Desktop/wgbs2.0/snakemake-test/genome/'
index = 'genome/'
rule all:
  input:
    expand("rawdata/{sample}_1.fq.gz",sample=SAMPLES),
    expand("rawdata/{sample}_2.fq.gz",sample=SAMPLES),
    expand("clean_fastq/{sample}_1.fq.gz",sample=SAMPLES),
    expand("clean_fastq/{sample}_2.fq.gz",sample=SAMPLES),
    expand("bam/{sample}_1_bismark_bt2_pe.bam",sample=SAMPLES),
    expand("bam/{sample}_1_bismark_bt2_PE_report.txt",sample=SAMPLES),
    expand("rmdup_bam/{sample}.deduplicated.bam",sample=SAMPLES),
    expand("rmdup_bam/{sample}.deduplication_report.txt",sample=SAMPLES),
    expand("cytosine_report/{sample}.deduplicated.CpG_report.txt.gz",sample=SAMPLES)
rule QC:
  input:
    raw_R1 = "rawdata/{sample}_1.fq.gz",
    raw_R2 = "rawdata/{sample}_2.fq.gz"
  output:
    clean_R1 = "clean_fastq/{sample}_1.fq.gz",
    clean_R2 = "clean_fastq/{sample}_2.fq.gz"
  threads: 4
  log:
    "clean_fastq/{sample}.html"
  shell:
    "fastp -w {threads} -i {input.raw_R1} -o {output.clean_R1} "
    "-I {input.raw_R2} -O {output.clean_R2} --detect_adapter_for_pe --html {log}"

rule Bismark_map:
  input:
    clean_R1 = "clean_fastq/{sample}_1.fq.gz",
    clean_R2 = "clean_fastq/{sample}_2.fq.gz"
  output:
    bam = "bam/{sample}_1_bismark_bt2_pe.bam",
    bam_report = "bam/{sample}_1_bismark_bt2_PE_report.txt"
  threads: 4
  log:
    "bam/{sample}_map_log.txt"
  shell:
    "bismark --parallel {threads} --non_directional --genome {index} -1 {input.clean_R1} -2 {input.clean_R2} -o ./bam/ 2>{log}"

rule samtools_sort:
  input:
    'bam/{sample}_1_bismark_bt2_pe.bam'
  output:
    'sortedbam/{sample}.bam'
  threads: 4
  shell:
    'samtools sort -n -@ {threads} -o {output} {input}'
    
rule Bismark_deduplicate:
  input:
    "sortedbam/{sample}.bam"
  output:
    rmdup_bam = "rmdup_bam/{sample}.deduplicated.bam",
    rmdup_report = "rmdup_bam/{sample}.deduplication_report.txt"
  log:
    "rmdup_bam/{sample}_rmdup_log.txt"
  params:
    output_filename = '{sample}',
    output_dir = 'rmdup_bam/',
  shell:
    "deduplicate_bismark -p {input} --output_dir {params.output_dir} -o {params.output_filename} 2>{log}"

rule Bismark_methylation_extractor:
  input:
    "rmdup_bam/{sample}.deduplicated.bam"
  output:
    "cytosine_report/{sample}.deduplicated.CpG_report.txt.gz"
  params:
    output_dir = 'cytosine_report/',
  threads: 4
  shell:
    "bismark_methylation_extractor --genome_folder {genome} --gzip --cytosine_report -p -o {params.output_dir} {input}"
