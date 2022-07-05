test_Peak <- paste0('test_Peak/',list.files('test_Peak/'))
test_peak_name <- stringr::str_extract(list.files('test_Peak/'),'.*(?=_rmdup_summits)')
all_ReproduciblePeaks <- indentifyReproduciblePeaks(summitFiles = test_Peak,summitNames = test_peak_name)
generate_gtf <- function(test_res){
  test_res <- as.data.frame(test_res)
  test2 <- data.frame(
    chr = test_res$seqnames,
    source = rep('.',nrow(test_res)),
    feature = rep('exon',nrow(test_res)),
    start = test_res$start,
    end = test_res$end,
    score = rep('.',nrow(test_res)),
    strand = rep('.',nrow(test_res)),
    frame = rep('.',nrow(test_res)),
    attribute = rep('gene_id',nrow(test_res)),
    name =  paste0(rep('Peak',nrow(test_res)),1:nrow(test_res))
  )
  return(test2)
}

all_ReproduciblePeaks_gtf <- generate_gtf(all_ReproduciblePeaks)
write.table(all_ReproduciblePeaks_gtf,file = 'all_ReproduciblePeaks.gtf',sep = "\t",col.names = F,row.names = F,quote = F)