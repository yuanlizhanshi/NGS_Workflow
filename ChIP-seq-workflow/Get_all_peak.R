library(tidyverse)
library(GenomicRanges)
generate_gtf <- function(peak_info) {
  starts <- peak_info$start + round(peak_info$length/2  - 251 )
  ends <-  peak_info$start + round(peak_info$length/2 + 249 )
  if (any(starts< 0)) {
    starts[which(starts< 0)] <- 1
    ends[which(starts< 0)] <- 501
  }
  test <- data.frame(
    chr = peak_info$chr,
    start =  starts,
    end = ends
  )
  test <- as.data.frame(GenomicRanges::reduce(makeGRangesFromDataFrame(test)))
  test2 <- data.frame(
    chr = test$seqnames,
    source = rep('.',nrow(test)),
    feature = rep('exon',nrow(test)),
    start = test$start,
    end = test$end,
    score = rep('.',nrow(test)),
    strand = rep('.',nrow(test)),
    frame = rep('.',nrow(test)),
    name = paste0('gene_id "',
                  paste0(rep('Peak',nrow(test)),1:nrow(test)),
                  '"')
  )
  return(test2)
}
sample <-  c("test1","test2")

all_peak <- map_dfr(paste0("peak/",sample,"/",paste0(sample,"_peaks.xls")),read.table,header = T)
write.table(all_peak,file = 'peak/all_peak_info.txt',sep = "\t",col.names = T,row.names = F,quote = F)
all_peak_gtf <- generate_gtf(all_peak)
write.table(all_peak_gtf,file = 'peak/all_peak_info.gtf',sep = "\t",col.names = F,row.names = F,quote = F)



