import glob
import os

directory_path = "libraries"

pattern = "*.csv"
library_files = []

for file_path in glob.glob(os.path.join(directory_path, pattern)):
    library_files.append(os.path.splitext(os.path.basename(file_path))[0])

cellranger = '~/software/cellranger-arc-2.0.2/bin/cellranger-arc'
cellranger_index = '/home/kyh/Desktop/hg38/CellrangerArc_ucsc_hg38'

# ~/software/cellranger-arc-2.0.2/bin/cellranger-arc count --id=R207 \
# --reference=/home/kyh/Desktop/hg38/CellrangerArc_ucsc_hg38 \
# --libraries=/data/gaowei/Islet_multiome_fastq/islet_multiome_out/libraries/R207_library.csv \
# --localcores=70 \
# --localmem=512


rule all:
	input:
		expand("libraries/{sample}.csv",sample=library_files),
		expand('{sample}/outs/summary.csv',sample=library_files)


rule cellranger:
	input:
		'libraries/{sample}.csv'
	threads: 70
	output:
		'{sample}/outs/summary.csv'
	shell:
		"""
		rm -rf {wildcards.sample}
		{cellranger} count --id={wildcards.sample} --reference={cellranger_index} --libraries={input} --localcores={threads} --localmem=512
		"""

    
