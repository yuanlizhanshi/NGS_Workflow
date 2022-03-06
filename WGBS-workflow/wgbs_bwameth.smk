SAMPLES = {'test_1','test_2'}
genome = 'genome/chr1_2.fa'

rule all:
  input:
    expand("rawdata/{sample}_1.fq.gz",sample=SAMPLES),
    expand("rawdata/{sample}_2.fq.gz",sample=SAMPLES),
    expand("clean_fastq/{sample}_1.fq.gz",sample=SAMPLES),
    expand("clean_fastq/{sample}_2.fq.gz",sample=SAMPLES),
    expand("rmdup_bam/{sample}_rmdup.bam",sample=SAMPLES),
    expand("rmdup_bam/{sample}_rmdup.bam.bai",sample=SAMPLES),
    expand("rmdup_bam/{sample}_rmdup.bam.bai",sample=SAMPLES),
    expand("methylKit/{sample}_CpG.methylKit",sample=SAMPLES),
    expand("cytosine_report/{sample}.cytosine_report.txt.gz",sample=SAMPLES)

rule QC:
  input:
    raw_R1 = "rawdata/{sample}_1.fq.gz",
    raw_R2 = "rawdata/{sample}_2.fq.gz"
  output:
    clean_R1 = "clean_fastq/{sample}_1.fq.gz",
    clean_R2 = "clean_fastq/{sample}_2.fq.gz"
  threads: 4
  shell:
    "fastp -w {threads} -i {input.raw_R1} -o {output.clean_R1} "
    "-I {input.raw_R2} -O {output.clean_R2}"

rule Bwameth_map:
  input:
    clean_R1 = "clean_fastq/{sample}_1.fq.gz",
    clean_R2 = "clean_fastq/{sample}_2.fq.gz"
  output:
    "sam/{sample}.sam"
  threads: 4
  shell:
    "bwameth.py --reference {genome} {input.clean_R1} {input.clean_R2} -t {threads} > {output}"

rule samtools_sort:
  input:
    'sam/{sample}.sam'
  output:
    'sortedbam/{sample}.bam'
  threads: 4
  shell:
    'samtools sort -@ {threads} -o {output} {input}'

rule samtools_remove_duplication:
  input:
    'sortedbam/{sample}.bam'
  output:
    bam = 'rmdup_bam/{sample}_rmdup.bam',
  shell:
    'samtools rmdup {input} {output.bam}'

rule samtools_index:
  input:
    'rmdup_bam/{sample}_rmdup.bam'   
  output:
    'rmdup_bam/{sample}_rmdup.bam.bai'
  shell:
    "samtools index {input}"

rule MethylDackel_extract_methylKit:
  input:
    'rmdup_bam/{sample}_rmdup.bam'
  output:
    "methylKit/{sample}_CpG.methylKit"
  params:
    output_prefix = 'methylKit/{sample}'
  shell:
    "MethylDackel extract -@ {threads} --methylKit {genome} {input} -o {params.output_prefix}"

rule MethylDackel_cytosine_report:
  input:
    'rmdup_bam/{sample}_rmdup.bam'
  output:
    "cytosine_report/{sample}.cytosine_report.txt"
  params:
    output_prefix = 'cytosine_report/{sample}'
  shell:
    "MethylDackel extract -@ {threads} --CHH --CHG --cytosine_report {genome} {input} -o {params.output_prefix}"

rule Bgzip_cytosine_report:
  input:
    "cytosine_report/{sample}.cytosine_report.txt"
  output:
    "cytosine_report/{sample}.cytosine_report.txt.gz"
  threads: 4
  shell:
    "bgzip -@ {threads} {input}"
