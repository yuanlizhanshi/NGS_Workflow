import os
import glob
import re

directory_path = "sra_file"
SAMPLES = []
files = glob.glob(directory_path + '/*')
for i in files:
    SAMPLES.append(os.path.splitext(os.path.basename(i))[0])

rule all:
    input:
        expand("sra_file/{sample}",sample=SAMPLES),
        expand("fastq/{sample}_1.fastq.gz",sample=SAMPLES)

rule sra2fastqs:
    input:
        file= 'sra_file/{sample}'
    output:
        file= "fastq/{sample}_1.fastq.gz"
    shell:
        """
        fastq-dump --split-files --gzip --defline-qual '+' --defline-seq '@$ac-$si/$ri' -O ./fastq/ {input}
        """

        