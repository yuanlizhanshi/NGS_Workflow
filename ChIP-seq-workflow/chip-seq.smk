SAMPLES = {'test1','test2'}
TYPES = {"input","IP"}
genome = 'genome/chr1_2.fa'
Bowtie2_index = "genome/chr1_2"

rule all:
  input:
    expand("rawdata/{sample}_{type}_1.fastq.gz",sample=SAMPLES,type = TYPES),
    expand("rawdata/{sample}_{type}_1.fastq.gz",sample=SAMPLES,type = TYPES),
    expand("clean_fastq/{sample}_{type}_1.fq.gz",sample=SAMPLES,type = TYPES),
    expand("clean_fastq/{sample}_{type}_2.fq.gz",sample=SAMPLES,type = TYPES),
    expand("rmdup_bam/{sample}_{type}_rmdup.bam",sample=SAMPLES,type = TYPES),
    expand("peak/{sample}/{sample}_peaks.narrowPeak",sample=SAMPLES),
    expand("peak/{sample}/{sample}_peaks.xls",sample=SAMPLES),
    expand("peak/{sample}/{sample}_summits.bed",sample=SAMPLES)

rule QC:
  input:
    raw_R1 = "rawdata/{sample}_{type}_1.fastq.gz",
    raw_R2 = "rawdata/{sample}_{type}_2.fastq.gz"
  output:
    clean_R1 = "clean_fastq/{sample}_{type}_1.fq.gz",
    clean_R2 = "clean_fastq/{sample}_{type}_2.fq.gz"
  threads: 4
  shell:
    "fastp -w {threads} -i {input.raw_R1} -o {output.clean_R1} "
    "-I {input.raw_R2} -O {output.clean_R2} --detect_adapter_for_pe"

rule Bowtie2_map:
  input:
    clean_R1 = "clean_fastq/{sample}_{type}_1.fq.gz",
    clean_R2 = "clean_fastq/{sample}_{type}_2.fq.gz"
  output:
    "sam/{sample}_{type}.sam"
  threads: 4
  log:
    "sam/{sample}_{type}_mapping_log.txt"
  shell:
    "bowtie2 -p {threads} --local -N 1 -1 {input.clean_R1} -2 {input.clean_R2} -x {Bowtie2_index} -S {output} 2>{log}"

rule samtools_sort:
  input:
    'sam/{sample}_{type}.sam'
  output:
    'sortedbam/{sample}_{type}.bam'
  threads: 4
  shell:
    'samtools sort -@ {threads} -o {output} {input}'

rule samtools_remove_duplication:
  input:
    'sortedbam/{sample}_{type}.bam'
  output:
    bam = 'rmdup_bam/{sample}_{type}_rmdup.bam',
  shell:
    'samtools rmdup {input} {output.bam}'

rule bam_index:
  input:
    'rmdup_bam/{sample}_{type}_rmdup.bam'
  output:
    'rmdup_bam/{sample}_{type}_rmdup.bam..bai'
  shell:
    "samtools index {input}"


rule bam_to_bigwig:
  input:
    bam = 'rmdup_bam/{sample}_{type}_rmdup.bam',
    bam_index = 'rmdup_bam/{sample}_{type}_rmdup.bam..bai'
  output:
    'rmdup_bam/{sample}_{type}_rmdup.bigwig'
  log:
    'bigwig/{sample}_{type}.bigwig.log'
  shell:
    'bamCoverage -bs 50 -b {input} -o {output} --normalizeUsing CPM 2>{log}'

rule Macs2_Peak_calling:
  input:
    IP_sample = "rmdup_bam/{sample}_IP_rmdup.bam",
    input_sample = "rmdup_bam/{sample}_input_rmdup.bam"
  output:
    Narrow_peak = "peak/{sample}/{sample}_peaks.narrowPeak",
    Peak= "peak/{sample}/{sample}_peaks.xls",
    Summits = "peak/{sample}/{sample}_summits.bed"
  params:
    output_prefix = "{sample}",
    output_dir = "./peak/{sample}/"
  shell:
    "macs2 callpeak -t {input.IP_sample} -c {input.input_sample} -f BAMPE -g 4.5e8 -n {params.output_prefix} -B -q 0.001 --outdir {params.output_dir}"

