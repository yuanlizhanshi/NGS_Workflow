import glob
import os

directory_path = "library"

pattern = "*.sh"
library_files = []

for file_path in glob.glob(os.path.join(directory_path, pattern)):
    library_files.append(os.path.splitext(os.path.basename(file_path))[0])


rule all:
	input:
		expand("library/{sample}.sh",sample=library_files),
		expand('{sample}/outs/summary.csv',sample=library_files)

rule cellranger:
	input:
		'library/{sample}.sh'
	threads: 10
	output:
		'{sample}/outs/summary.csv'
	shell:
		"""
		rm -rf {wildcards.sample}
		bash {input}
		"""

    
