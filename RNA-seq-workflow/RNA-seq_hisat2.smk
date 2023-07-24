SAMPLES = {'js','L10','Pps'}
genome = 'genome/chr1.fa'
gtf = 'genome/chr1.gtf'
index = 'genome/chr1_index'
rule all:
  input:
    expand("rawdata/{sample}_1.fq.gz",sample=SAMPLES),
    expand("rawdata/{sample}_2.fq.gz",sample=SAMPLES),
    expand("clean_fastq/{sample}_1.fq.gz",sample=SAMPLES),
    expand("clean_fastq/{sample}_2.fq.gz",sample=SAMPLES),
    expand("sortedbam/{sample}.bam",sample=SAMPLES),
    "counts.txt"

rule QC:
  input:
    raw_R1 = "rawdata/{sample}_1.fq.gz",
    raw_R2 = "rawdata/{sample}_2.fq.gz"
  output:
    clean_R1 = "clean_fastq/{sample}_1.fq.gz",
    clean_R2 = "clean_fastq/{sample}_2.fq.gz"
  log:
    "clean_fastq/{sample}.html"
  threads: 4
  shell:
    "fastp -w {threads} -i {input.raw_R1} -o {output.clean_R1} "
    "-I {input.raw_R2} -O {output.clean_R2} --detect_adapter_for_pe --html {log}"

rule hisat2_map:
  input:
    clean_R1 = "clean_fastq/{sample}_1.fq.gz",
    clean_R2 = "clean_fastq/{sample}_2.fq.gz"
  output:
    temp('sam/{sample}.sam')
  log:
    "sam/{sample}_mapping_log.txt"
  threads: 4
  shell:
    "hisat2 -p {threads} -x {index} --dta -1 {input.clean_R1} -2 {input.clean_R2} -S {output} 2>{log}"

rule samtools_sort:
  input:
    'sam/{sample}.sam'
  output:
    'sortedbam/{sample}.bam'
  threads: 4
  shell:
    'samtools sort -@ {threads} -o {output} {input}'

rule counts:
  input:
    gtf = {gtf},
    bam = expand('sortedbam/{sample}.bam',sample=SAMPLES)
  output:
    "counts.txt"
  threads: 4
  shell:
    "featureCounts -a {input.gtf} -o {output} -T {threads} {input.bam}"

#For paired reads
#rule counts:
#   input:
#     gtf = {gtf},
#     bam = expand('sortedbam/{sample}.bam',sample=SAMPLES)
#   output:
#     "counts.txt"
#   threads: 4
#   shell:
#     "featureCounts -p --countReadPairs -a {input.gtf} -o {output} -T {threads} {input.bam}"
