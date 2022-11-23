samples = {'JS','L10','big','small'}

genome = './genome/silkworm.fa'
GATK = '/home/kyh/software/gatk/gatk'

fai = open('/home/kyh/Desktop/refindex/silkworm.fa.fai','r')
chr = []
for line in fai:
  line = line.split('\t')
  chr.append(line[0])
fai.close()

#wildcard_constraints:
#    samples = "[A-Z0-9]",
#    chr= ".{11}"
    
rule all:
  input:
    expand("/home/kyh/Desktop/a42/silkdb3.0/rawdata/{sample}_1.fq.gz",sample=samples),
    expand("/home/kyh/Desktop/a42/silkdb3.0/rawdata/{sample}_2.fq.gz",sample=samples),
    expand("/home/kyh/Desktop/a42/silkdb3.0/clean_fastq/{sample}_1.fq.gz",sample=samples),
    expand("/home/kyh/Desktop/a42/silkdb3.0/clean_fastq/{sample}_2.fq.gz",sample=samples),
    expand("markdup_bam/{sample}_markdup.bam",sample=samples),
    expand("markdup_bam/{sample}_markdup.bam.bai",sample=samples),
    expand('last_bam/{sample}_last.bam',sample= samples),
    expand('last_gvcf/{sample}_last.gvcf',sample= samples),
    expand("vcf/all_sample.vcf.gz")
rule QC:
  input:
    raw_R1 = "/home/kyh/Desktop/a42/silkdb3.0/rawdata/{sample}_1.fq.gz",
    raw_R2 = "/home/kyh/Desktop/a42/silkdb3.0/rawdata/{sample}_2.fq.gz"
  output:
    clean_R1 = "/home/kyh/Desktop/a42/silkdb3.0/clean_fastq/{sample}_1.fq.gz",
    clean_R2 = "/home/kyh/Desktop/a42/silkdb3.0/clean_fastq/{sample}_2.fq.gz"
  threads: 16
  log:
    "/home/kyh/Desktop/a42/silkdb3.0/clean_fastq/{sample}.html"   
  shell:
    "fastp -w {threads} -i {input.raw_R1} -o {output.clean_R1} "
    "-I {input.raw_R2} -O {output.clean_R2} --detect_adapter_for_pe --html {log}"
rule Bwa_map:
  input:
    clean_R1 = "/home/kyh/Desktop/a42/silkdb3.0/clean_fastq/{sample}_1.fq.gz",
    clean_R2 = "/home/kyh/Desktop/a42/silkdb3.0/clean_fastq/{sample}_2.fq.gz"
  output:
    temp("sam/{sample}.sam")
  threads: 20
  shell:
    "bwa mem -t {threads} -R '@RG\\tID:{wildcards.sample}\\tPL:illuminar\\tSM:{wildcards.sample}' "
    "{genome} {input.clean_R1} {input.clean_R2} > {output}"

rule samtools_sort:
  input:
    'sam/{sample}.sam'
  output:
    temp('sortedbam/{sample}.bam')
  threads: 20
  shell:
    'samtools sort -@ {threads} -o {output} {input}'

rule GATK_MarkDuplicates:
  input:
    'sortedbam/{sample}.bam'
  output:
    bam = 'markdup_bam/{sample}_markdup.bam',
    metrix = temp('markdup_bam/{sample}.metrix')
  shell:
    '{GATK} --java-options "-Xmx100G" MarkDuplicates -I {input} -M {output.metrix} -O {output.bam}'

rule samtools_index:
  input:
    'markdup_bam/{sample}_markdup.bam'
  output:
    'markdup_bam/{sample}_markdup.bam.bai'
  shell:
    "samtools index {input}"
    
rule Variant_calling_by_chromosomes:
    input:
        bam = 'markdup_bam/{sample}_markdup.bam',
        bai = 'markdup_bam/{sample}_markdup.bam.bai'
    output:
        bam = temp('temp_bam/{sample}.{chromosome}.bam'),
        gvcf = temp('temp_gvcf/{sample}.{chromosome}.gvcf')
    shell:
        """
        {GATK} --java-options "-Xmx100G -XX:ParallelGCThreads=4" HaplotypeCaller -R {genome} --emit-ref-confidence GVCF -I {input.bam} -O {output.gvcf} -bamout {output.bam} -L {wildcards.chr}
        """
rule merge_bam_by_chromosomes:
    input:
      bam = expand('temp_bam/{{sample}}.{chromosome}.bam', sample= samples, chromosome= chr)
    output:
      'last_bam/{sample}_last.bam'
    params:
      extra = lambda wildcards, input: ' '.join(input.bam)
    shell:
      "samtools merge {output} {params.extra}"

rule merge_gvcf_by_chromosomes:
    input:
      gvcf = expand('temp_gvcf/{{sample}}.{chromosome}.gvcf', sample= samples, chromosome= chr)
    output:
      'last_gvcf/{sample}_last.gvcf'
    params:
      extra = lambda wildcards, input: ' -V '.join(input.gvcf)
    shell:
      "{GATK} CombineGVCFs -R {genome} -V {params.extra} -O {output}"

rule Combine_gvcf:
  input:
    gvcf = expand("last_gvcf/{sample}_last.gvcf",sample = samples)
  output:
    temp("gvcf/all_sample.gvcf.gz")
  params:
    extra = lambda wildcards, input: ' -V '.join(input.gvcf)
  shell:
    "{GATK} CombineGVCFs -R {genome} -V {params.extra} -O {output}"

rule Genotype_gvcf:
  input:
    "gvcf/all_sample.gvcf.gz"
  output:
    "vcf/all_sample.vcf.gz"
  shell:
    '{GATK} --java-options "-Xmx100G" GenotypeGVCFs -R {genome} -V {input} -O {output}'


rule Extract_SNP:
  input:
    "vcf/all_sample.vcf.gz"
  output:
    "vcf/all_sample_snp.vcf.gz"
  shell:
    '{GATK} --java-options "-Xmx100G" SelectVariants -R {genome} '
    '-V {input} -O {output} --select-type-to-include SNP'
    
rule Extract_indel:
  input:
    "vcf/all_sample.vcf.gz"
  output:
    "vcf/all_sample_indel.vcf.gz"
  shell:
    '{GATK} --java-options "-Xmx100G" SelectVariants -R {genome} '
    '-V {input} -O {output} --select-type-to-include INDEL'
    
