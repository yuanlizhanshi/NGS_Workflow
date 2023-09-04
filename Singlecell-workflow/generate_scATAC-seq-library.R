library(tidyverse)

libraries <- list.files('./hpapdata2/',pattern = 'HPAP')
generate_library <- function(x) {
  fastqs <- paste0('/data/kyh/hpap/hpapdata2/',x)
  sample <- unique(str_extract(list.files(paste0('hpapdata2/',x)),'.*(?=_S)'))
  library_type <- rep('Chromatin Accessibility',length(sample))
  df <- tibble(
    fastqs = fastqs,
    sample = sample,
    library_type = library_type
  )
  df <- distinct_all(df)
  return(df)
}

hpap_all <- map_dfr(libraries,generate_library) 
data.table::fwrite(hpap_all,file = 'hpapdata2/hpap_sample_info.txt',sep = '\t')

library(glue)
cellranger = '~/software/cellranger-atac-2.1.0/cellranger-atac'
cellranger_index = '/home/kyh/Desktop/hg38/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/'
for (i in 1:nrow(hpap_all)) {
    samplename = hpap_all$sample[i]
    fastq = hpap_all$fastqs[i]
    command = glue("{cellranger} count --id={samplename} --reference={cellranger_index} --fastqs={fastq} --sample={samplename} --localcores=10 --localmem=512")
    write_lines(command,file = glue('hpap_out/library/{samplename}.sh'))
}

cellranger_aggr <- tibble(
    library_id = hpap_all$sample,
    fragments = paste0('/data/kyh/hpap/hpap_out/',hpap_all$sample,'/outs/fragments.tsv.gz'),
    cells = paste0('/data/kyh/hpap/hpap_out/',hpap_all$sample,'/outs/singlecell.csv')
)

write_csv(cellranger_aggr,'/data/kyh/hpap/hpap_out/cellrangera_hpap_atac_aggr_libraries.csv')

# ~/software/cellranger-atac-2.1.0/cellranger-atac aggr \
# --id=hpap_atac_aggr \
# --csv=cellrangera_hpap_atac_aggr_libraries.csv \
# --normalize=depth \
# --localcores=40 \
# --reference=/home/kyh/Desktop/hg38/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/