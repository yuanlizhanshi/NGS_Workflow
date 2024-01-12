from pathlib import Path
import re
import pandas as pd
all_file = Path("rawdata/").glob('*fastq.gz')
all_file_list = [i for i in all_file]

SAMPLES = []
for i in all_file_list:
    file = re.search('.*(?=_read\\d.fastq)',i.stem).group()
    if file not in SAMPLES:
        SAMPLES.append(file)

hg38 = "~/Desktop/hg38/hg38_bowtie2/hg38.fa"
Bowtie2_index = "~/Desktop/hg38/hg38_bowtie2/hg38"
gtf = '~/Desktop/hg38/hg38_transcript/genes.gtf'
genome_black_list = "/data/kyh/NB4_GRN/basic_GRN/NB4_ABC/reference/hg38_blacklist.bed"
r_path = "~/miniconda3/envs/r43/bin/Rscript"
# my_dict = {}
# def get_mapping_rate(file):
#     with open(file,'r') as f:
#         for line in f:
#             if line.find('overall') != -1:
#                 res = re.search('.*(?=\\so)',line).group()
#                 sample = re.search('.*(?=_mapping_log)',file.stem).group()
#                 my_dict[sample] = res
#             
# for i in all_file_list:
#     get_mapping_rate(i)
# df = pd.DataFrame(my_dict.items(),columns=['Sample','Mapping rate'])
# df.to_csv('sam/all_sample_mapping_rate.txt',sep='\t',index = False)

rule all:
    input:
        expand("rawdata/{sample}_read1.fastq.gz",sample=SAMPLES),
        expand("rawdata/{sample}_read2.fastq.gz",sample=SAMPLES),
        expand("clean_fastq/{sample}_1.fq.gz",sample=SAMPLES),
        expand("clean_fastq/{sample}_2.fq.gz",sample=SAMPLES),
        expand("rmdup_bam/{sample}_rmdup.bam",sample=SAMPLES),
        expand("rmdup_bam/{sample}_rmdup.bam.bai",sample=SAMPLES),
        expand("fragment_distuibution/{sample}_rmdup.pdf",sample=SAMPLES),
        expand("peak/{sample}/{sample}_peaks.narrowPeak",sample=SAMPLES),
        expand("peak/{sample}/{sample}_peaks.xls",sample=SAMPLES),
        expand("peak/{sample}/{sample}_summits.bed",sample=SAMPLES),
        expand("peak/{sample}/{sample}_peaks_finalhits.txt",sample=SAMPLES),
        expand("TOBIAS/{sample}_ATACorrect/{sample}_rmdup_corrected.bw",sample=SAMPLES),
        expand("TOBIAS/{sample}_ATACorrect/{sample}_footprints.bw",sample=SAMPLES)
        
rule QC:
    input:
        raw_R1 = "rawdata/{sample}_read1.fastq.gz",
        raw_R2 = "rawdata/{sample}_read2.fastq.gz"
    output:
        clean_R1 = "clean_fastq/{sample}_1.fq.gz",
        clean_R2 = "clean_fastq/{sample}_2.fq.gz"
    threads: 10
    log:
        "clean_fastq/{sample}.html"
    shell:
        "fastp -w {threads} -i {input.raw_R1} -o {output.clean_R1} "
        "-I {input.raw_R2} -O {output.clean_R2} --detect_adapter_for_pe --html {log} "

    
rule Bowtie2_map:
    input:
        clean_R1 = "clean_fastq/{sample}_1.fq.gz",
        clean_R2 = "clean_fastq/{sample}_2.fq.gz"
    output:
        temp("sam/{sample}.sam")
    threads: 30
    log:
        "sam/{sample}_mapping_log.txt"
    shell:
        "bowtie2 -p {threads} -X 2000 -1 {input.clean_R1} -2 {input.clean_R2} -x {Bowtie2_index} -S {output} 2>{log}"

rule samtools_sort:
    input:
        'sam/{sample}.sam'
    output:
        'sortedbam/{sample}.bam'
    threads: 30
    shell:
        'samtools sort -@ {threads} -o {output} {input}'

rule samtools_remove_duplication:
    input:
        'sortedbam/{sample}.bam'
    output:
        bam = 'rmdup_bam/{sample}_rmdup.bam',
    shell:
        'samtools rmdup {input} {output.bam}'

rule bam_index:
    input:
        'rmdup_bam/{sample}_rmdup.bam'
    output:
        'rmdup_bam/{sample}_rmdup.bam.bai'
    shell:
        "samtools index {input}"

rule bam_qc:
    input:
        'rmdup_bam/{sample}_rmdup.bam'
    output:
        'fragment_distuibution/{sample}_rmdup.pdf'
    shell:
        "{r_path} atac_qc.R --bam {input}"
        
  
rule Macs2_Peak_calling:
    input:
        "rmdup_bam/{sample}_rmdup.bam"
    output:
        Narrow_peak = "peak/{sample}/{sample}_peaks.narrowPeak",
        Peak= "peak/{sample}/{sample}_peaks.xls",
        Summits = "peak/{sample}/{sample}_summits.bed"
    params:
        output_prefix = "{sample}",
        output_dir = "./peak/{sample}/"
    shell:
        "macs2 callpeak -t {input} -f BAMPE --nomodel --shift -75 --extsize 150 -g hs -n {params.output_prefix} -B -q 0.01 --outdir {params.output_dir}"

rule peak_annotation:
    input:
        "peak/{sample}/{sample}_peaks.narrowPeak"
    output:
        "peak/{sample}/{sample}_peaks_finalhits.txt",
    params:
        output_dir = "./peak/{sample}/"
    shell:
        "uropa --bed {input} --gtf {gtf} --show_attributes gene_id gene_name --feature_anchor start --distance 3000 3000 --feature gene -o {params.output_dir}"

rule ATACorrect:
    input:
        bam =  "rmdup_bam/{sample}_rmdup.bam",
        peak = "peak/{sample}/{sample}_peaks.narrowPeak"
    output:
        "TOBIAS/{sample}_ATACorrect/{sample}_rmdup_corrected.bw",
    threads: 40
    params:
        output_dir = "TOBIAS/{sample}_ATACorrect"
    shell:
        '''
        TOBIAS ATACorrect --bam {input.bam} --genome {hg38} --peaks {input.peak} --blacklist {genome_black_list} --outdir {params.output_dir} --cores {threads}
        '''

rule FootprintScores:
    input:
        bigwig =  "TOBIAS/{sample}_ATACorrect/{sample}_rmdup_corrected.bw",
        peak = "peak/{sample}/{sample}_peaks.narrowPeak"
    output:
        "TOBIAS/{sample}_ATACorrect/{sample}_footprints.bw",
    threads: 40
    shell:
        '''
        TOBIAS FootprintScores --signal {input.bigwig} --regions {input.peak} --output {output} --cores {threads}
        '''