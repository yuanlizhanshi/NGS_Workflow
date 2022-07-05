SAMPLES = {'qinJS_test','qinL10_test'}
genome = './genome/chr1_2.fa'
GATK = '/home/kyh/Desktop/gatk/gatk'

rule all:
  input:
    expand("rawdata/{sample}_1.fastq.gz",sample=SAMPLES),
    expand("rawdata/{sample}_2.fastq.gz",sample=SAMPLES),
    expand("clean_fastq/{sample}_1.fq.gz",sample=SAMPLES),
    expand("clean_fastq/{sample}_2.fq.gz",sample=SAMPLES),
    expand("markdup_bam/{sample}_markdup.bam",sample=SAMPLES),
    expand("markdup_bam/{sample}_markdup.bam.bai",sample=SAMPLES),
    expand("vcf/all_sample.vcf.gz")
rule QC:
  input:
    raw_R1 = "rawdata/{sample}_1.fastq.gz",
    raw_R2 = "rawdata/{sample}_2.fastq.gz"
  output:
    clean_R1 = "clean_fastq/{sample}_1.fq.gz",
    clean_R2 = "clean_fastq/{sample}_2.fq.gz"
  threads: 4
  log:
    "clean_fastq/{sample}.html"   
  shell:
    "fastp -w {threads} -i {input.raw_R1} -o {output.clean_R1} "
    "-I {input.raw_R2} -O {output.clean_R2} --detect_adapter_for_pe --html {log}"
rule Bwa_map:
  input:
    clean_R1 = "clean_fastq/{sample}_1.fq.gz",
    clean_R2 = "clean_fastq/{sample}_2.fq.gz"
  output:
    "sam/{sample}.sam"
  threads: 4
  shell:
    "bwa mem -t {threads} -R '@RG\\tID:{wildcards.sample}\\tPL:illuminar\\tSM:{wildcards.sample}' "
    "{genome} {input.clean_R1} {input.clean_R2} > {output}"

rule samtools_sort:
  input:
    'sam/{sample}.sam'
  output:
    'sortedbam/{sample}.bam'
  threads: 4
  shell:
    'samtools sort -@ {threads} -o {output} {input}'

rule GATK_MarkDuplicates:
  input:
    'sortedbam/{sample}.bam'
  output:
    bam = 'markdup_bam/{sample}_markdup.bam',
    metrix = 'markdup_bam/{sample}.metrix'
  shell:
    '{GATK} MarkDuplicates -I {input} -M {output.metrix} -O {output.bam}'

rule samtools_index:
  input:
    'markdup_bam/{sample}_markdup.bam'
  output:
    'markdup_bam/{sample}_markdup.bam.bai'
  shell:
    "samtools index {input}"

rule Variant_calling:
  input:
    'markdup_bam/{sample}_markdup.bam'
  output:
    'gvcf/{sample}.g.vcf.gz',
  shell:
    '{GATK} HaplotypeCaller -R {genome} --emit-ref-confidence GVCF -I {input} -O {output}'

rule Combine_gvcf:
  input:
    gvcf = expand("gvcf/{sample}.g.vcf.gz",sample=SAMPLES)
  output:
    "vcf/all_sample.gvcf.gz"
  params:
    extra = lambda wildcards, input: ' -V '.join(input.gvcf)
  shell:
    "{GATK} CombineGVCFs -R {genome} -V {params.extra} -O {output}"

rule Genotype_gvcf:
  input:
    "vcf/all_sample.gvcf.gz"
  output:
    "vcf/all_sample.vcf.gz"
  shell:
    "{GATK} GenotypeGVCFs -R {genome} -V {input} -O {output}"
