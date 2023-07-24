SAMPLES = {'js','L10','Pps'}
STAR_index = '~/refindex/STARindex'
genome = 'genome/chr1.fa'
gtf = 'genome/chr1.gtf'
rule all:
  input:
    expand("raw_fastq/{sample}_1.fq.gz",sample=SAMPLES),
    expand("raw_fastq/{sample}_2.fq.gz",sample=SAMPLES),
    expand("clean_fastq/{sample}_1.fq.gz",sample=SAMPLES),
    expand("clean_fastq/{sample}_2.fq.gz",sample=SAMPLES),
    expand("bam/{sample}/Aligned.sortedByCoord.out.bam",sample=SAMPLES),
    'counts.txt'
rule QC:
  input:
    raw_R1 = "raw_fastq/{sample}_1.fq.gz",
    raw_R2 = "raw_fastq/{sample}_2.fq.gz"
  output:
    clean_R1 = "clean_fastq/{sample}_1.fq.gz",
    clean_R2 = "clean_fastq/{sample}_2.fq.gz"
  threads: 4
  log:
    "clean_fastq/{sample}.html"
  shell:
    "fastp -w {threads} -i {input.raw_R1} -o {output.clean_R1} "
    "-I {input.raw_R2} -O {output.clean_R2} --detect_adapter_for_pe --html --html {log}"

rule STAR_map:
  input:
    clean_R1 = "clean_fastq/{sample}_1.fq.gz",
    clean_R2 = "clean_fastq/{sample}_2.fq.gz"
  output:
    "bam/{sample}/Aligned.sortedByCoord.out.bam"
  threads: 4
  shell:
    "STAR --genomeDir {STAR_index} --twopassMode Basic "
    "--runThreadN {threads} --readFilesCommand zcat "
    "--readFilesIn {input.clean_R1} {input.clean_R2} "
    "--outSAMtype BAM SortedByCoordinate --sjdbOverhang 149 "
    "--outFilterIntronMotifs RemoveNoncanonical "
    "--outFileNamePrefix ./bam/{wildcards.sample}/"
rule rename:
  input:
    bam = "bam/{sample}/Aligned.sortedByCoord.out.bam"
  output:
    "bam/{sample}.bam"
  shell:
    "mv {input.bam} bam/{wildcards.sample}.bam;"

rule counts:
  input:
    gtf = {gtf},
    bam = expand('bam/{sample}.bam',sample=SAMPLES)
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
