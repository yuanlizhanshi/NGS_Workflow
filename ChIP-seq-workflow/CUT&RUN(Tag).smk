from pathlib import Path

Bowtie2_index = "~/Desktop/hg38/hg38_bowtie2/hg38"

all_file = Path("rawdata/").glob('*clean.fq.gz')
all_file_list = [i for i in all_file]

SAMPLES = []
for i in all_file_list:
    file = re.search('.*(?=_\\d.clean)',i.stem).group()
    if file not in SAMPLES:
        SAMPLES.append(file)
print(SAMPLES)

rule all:
    input:
        expand("rawdata/{sample}_1.clean.fq.gz",sample=SAMPLES),
        expand("rawdata/{sample}_2.clean.fq.gz",sample=SAMPLES),
        expand("clean_fastq/{sample}_1.fq.gz",sample=SAMPLES),
        expand("clean_fastq/{sample}_2.fq.gz",sample=SAMPLES),
        expand("rmdup_bam/{sample}_chrM_removed_dupRemoved.bam",sample=SAMPLES),
        expand("rmdup_bam/{sample}_chrM_removed_dupRemoved.bam.bai",sample=SAMPLES),
        expand("bigwig/{sample}.bigwig",sample=SAMPLES),
        expand("narrow_peak/{sample}/{sample}_peaks.xls",sample=SAMPLES),
        expand("broad_peak/{sample}/{sample}_peaks.xls",sample=SAMPLES)

rule QC:
    input:
        raw_R1 = "rawdata/{sample}_1.clean.fq.gz",
        raw_R2 = "rawdata/{sample}_2.clean.fq.gz"
    output:
        clean_R1 = "clean_fastq/{sample}_1.fq.gz",
        clean_R2 = "clean_fastq/{sample}_2.fq.gz"
    threads: 10
    log:
        "clean_fastq/{sample}.html"
    shell:
        """
        fastp -w {threads} -i {input.raw_R1} -o {output.clean_R1} -I {input.raw_R2} -O {output.clean_R2} --detect_adapter_for_pe --html {log}
        rm fastp.json
        """

    
rule Bowtie2_map:
    input:
        clean_R1 = "clean_fastq/{sample}_1.fq.gz",
        clean_R2 = "clean_fastq/{sample}_2.fq.gz"
    output:
        temp("sam/{sample}.sam")
    threads: 40
    log:
        "sam/{sample}_mapping_log.txt"
    shell:
        "bowtie2 -p {threads} --very-sensitive -X 2000 -1 {input.clean_R1} -2 {input.clean_R2} -x {Bowtie2_index} -S {output} 2>{log}"

rule samtools_sort:
    input:
        'sam/{sample}.sam'
    output:
        temp('bam/{sample}_sorted.bam')
    threads: 40
    shell:
        'samtools sort -@ {threads} -o {output} {input}'

rule samtools_remove_chrM:
    input:
        'bam/{sample}_sorted.bam'
    output:
        temp('bam/{sample}_sorted_chrM_removed.bam')
    threads: 40
    shell:
        'samtools view -h -q 30 -F 1804 -f 2 {input} |grep -v chrM | samtools sort -@ {threads} -O bam -o {output}'


rule samtools_remove_duplication:
    input:
        'bam/{sample}_sorted_chrM_removed.bam'
    output:
        'rmdup_bam/{sample}_chrM_removed_dupRemoved.bam',
    shell:
        'samtools rmdup {input} {output}'



rule bam_index:
    input:
        'rmdup_bam/{sample}_chrM_removed_dupRemoved.bam'
    output:
        'rmdup_bam/{sample}_chrM_removed_dupRemoved.bam.bai'
    shell:
        "samtools index {input}"


rule bam_to_bigwig:
    input:
        bam = 'rmdup_bam/{sample}_chrM_removed_dupRemoved.bam',
        bam_index = 'rmdup_bam/{sample}_chrM_removed_dupRemoved.bam.bai'
    output:
        'bigwig/{sample}.bigwig'
    threads: 20
    log:
        'bigwig/{sample}.bigwig.log'
    shell:
        'bamCoverage -bs 10 -p {threads} -b {input.bam} -o {output} --normalizeUsing RPGC --effectiveGenomeSize 2652783500 --ignoreForNormalization chrX 2>{log}'
       
rule Macs2_Narrow_Peak_calling:
    input:
        "rmdup_bam/{sample}_chrM_removed_dupRemoved.bam"
    output:
        Peak= "narrow_peak/{sample}/{sample}_peaks.xls",
    params:
        output_prefix = "{sample}",
        output_dir = "./narrow_peak/{sample}/"
    shell:
        "macs2 callpeak -t {input} -f BAMPE -n {params.output_prefix} -g mm --keep-dup 1 --cutoff-analysis --outdir {params.output_dir}"
        
rule Macs2_broad_Peak_calling:
    input:
        "rmdup_bam/{sample}_chrM_removed_dupRemoved.bam"
    output:
        Peak= "broad_peak/{sample}/{sample}_peaks.xls"
    params:
        output_prefix = "{sample}",
        output_dir = "./broad_peak/{sample}/"
    shell:
        "macs2 callpeak -t {input} -f BAMPE -n {params.output_prefix} -g mm --keep-dup 1 --cutoff-analysis --broad --outdir {params.output_dir}"      
        