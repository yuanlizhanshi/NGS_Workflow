import glob
import os

directory_path = "sra"


sra_files = []
files = glob.glob(directory_path + '/*')
for i in files:
    sra_files.append(os.path.basename(i))

rule all:
  input:
    expand("sra/{sample}",sample=sra_files),
    expand("fastq/{sample}_1.fastq.gz",sample=sra_files)
    
rule fastq_dump:
    input:
        "sra/{sample}"
    output:
        'fastq/{sample}_1.fastq.gz'
    threads: 100
    shell:
        "fastq-dump --gzip --split-files --defline-qual '+' --defline-seq '@$ac-$si/$ri' -O ./fastq/ {input}"
