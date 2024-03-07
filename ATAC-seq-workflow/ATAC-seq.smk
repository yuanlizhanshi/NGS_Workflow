SAMPLES = {"test1","test2"}

genome = 'genome/chr1_2.fa'
Bowtie2_index = "genome/chr1_2"
gtf = 'genome.gtf'

rule all:
	input:
		expand("rawdata/{sample}_1.fastq.gz",sample=SAMPLES),
		expand("rawdata/{sample}_2.fastq.gz",sample=SAMPLES),
		expand("clean_fastq/{sample}_1.fq.gz",sample=SAMPLES),
		expand("clean_fastq/{sample}_2.fq.gz",sample=SAMPLES),
		expand("rmdup_bam/{sample}_rmdup.bam",sample=SAMPLES),
		expand("peak/{sample}/{sample}_peaks.narrowPeak",sample=SAMPLES),
		expand("peak/{sample}/{sample}_peaks.xls",sample=SAMPLES),
		expand("peak/{sample}/{sample}_summits.bed",sample=SAMPLES),
		expand("peak/{sample}/{sample}_peaks_finalhits.txt",sample=SAMPLES)
    

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
		"-I {input.raw_R2} -O {output.clean_R2} --detect_adapter_for_pe --html {log} "

    
rule Bowtie2_map:
	input:
		clean_R1 = "clean_fastq/{sample}_1.fq.gz",
		clean_R2 = "clean_fastq/{sample}_2.fq.gz"
	output:
		temp("sam/{sample}.sam")
	threads: 4
	log:
		"sam/{sample}_mapping_log.txt"
	shell:
		"bowtie2 -p {threads} -N 1 -X 2000 -1 {input.clean_R1} -2 {input.clean_R2} -x {Bowtie2_index} -S {output} 2>{log}"

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

rule bam_index:
	input:
		'rmdup_bam/{sample}_rmdup.bam'
	output:
		'rmdup_bam/{sample}_rmdup.bam.bai'
	shell:
		"samtools index {input}"

rule bam_to_bigwig:
	input:
		bam = 'rmdup_bam/{sample}_rmdup.bam',
		bam_index = 'rmdup_bam/{sample}_rmdup.bam.bai'
	output:
		'bigwig/{sample}.bigwig'
	threads: 20
	log:
		'bigwig/{sample}.bigwig.log'
	shell:
		'bamCoverage -bs 20 -p {threads} -b {input.bam} -o {output} --normalizeUsing CPM 2>{log}'
       

rule Macs3_Peak_calling:
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
		"macs3 hmmratac -i {input} -n {params.output_prefix} --outdir {params.output_dir}"

rule peak_annotation:
    input:
        "peak/{sample}/{sample}_peaks.narrowPeak"
    output:
        "peak/{sample}/{sample}_peaks_finalhits.txt",
    params:
        output_dir = "./peak/{sample}/"
    shell:
        "uropa --bed {input} --gtf {gtf} --show_attributes gene_id gene_name --feature_anchor start --distance 3000 3000 --feature gene -o {params.output_dir}"