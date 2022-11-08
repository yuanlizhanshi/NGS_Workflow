library(tidyverse)
library(GenomicRanges)
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
indentifyReproduciblePeaks <- function (summitFiles = NULL, summitNames = NULL, reproducibility = 0.51, 
                                        extendSummits = 250, prefix = NULL) 
{
  .getQuantiles <- function (v = NULL, len = length(v)) 
  {
    if (length(v) < len) {
      v2 <- rep(0, len)
      v2[seq_along(v)] <- v
    }
    else {
      v2 <- v
    }
    p <- trunc(rank(v2))/length(v2)
    if (length(v) < len) {
      p <- p[seq_along(v)]
    }
    return(p)
  }
  nonOverlappingGR <- function (gr = NULL, by = "score", decreasing = TRUE) 
  {
    
    .clusterGRanges <- function(gr = NULL, filter = TRUE, by = "score", 
                                decreasing = TRUE) {
      gr <- sort(sortSeqlevels(gr))
      r <- GenomicRanges::reduce(gr, min.gapwidth = 0L, ignore.strand = TRUE)
      o <- findOverlaps(gr, r, ignore.strand = TRUE)
      mcols(gr)$cluster <- subjectHits(o)
      gr <- gr[order(mcols(gr)[, by], decreasing = decreasing), 
      ]
      gr <- gr[!duplicated(mcols(gr)$cluster), ]
      gr <- sort(sortSeqlevels(gr))
      mcols(gr)$cluster <- NULL
      return(gr)
    }
    i <- 0
    grConverge <- gr
    while (length(grConverge) > 0) {
      i <- i + 1
      grSelect <- .clusterGRanges(gr = grConverge, filter = TRUE, 
                                  by = by, decreasing = decreasing)
      grConverge <- subsetByOverlaps(grConverge, grSelect, 
                                     invert = TRUE, ignore.strand = TRUE)
      if (i == 1) {
        grAll <- grSelect
      }
      else {
        grAll <- c(grAll, grSelect)
      }
    }
    grAll <- sort(sortSeqlevels(grAll))
    return(grAll)
  }
  summit2gr <- function(summitfile){
    df <- read.table(summitfile)
    colnames(df)  <- c('seqnames','start','end','peak','score')
    gr <- makeGRangesFromDataFrame(df)
    gr$score <- df$score
    return(gr)
  }
  nonOverlapPassES <- tryCatch({
    summits <- lapply(seq_along(summitFiles), function(x) {
      grx <- summit2gr(summitFiles[x])
      grx$GroupReplicate <- paste0(summitNames[x])
      grx
    })
    summits <- Reduce("c", as(summits, "GRangesList"))
    
    extendedSummits <- resize(summits, extendSummits * 2 + 
                                1, "center")
    extendedSummits <- lapply(split(extendedSummits, extendedSummits$GroupReplicate), 
                              function(x) {
                                nonES <- nonOverlappingGR(x, by = "score", decreasing = TRUE)
                                nonES$replicateScoreQuantile <- round(.getQuantiles(nonES$score), 
                                                                      3)
                                nonES
                              })
    extendedSummits <- Reduce("c", as(extendedSummits, "GRangesList"))
    
    nonOverlapES <- nonOverlappingGR(extendedSummits, by = "replicateScoreQuantile", 
                                     decreasing = TRUE)
    
    overlapMat <- lapply(split(extendedSummits, extendedSummits$GroupReplicate), 
                         function(x) {
                           overlapsAny(nonOverlapES, x)
                         }) %>% Reduce("cbind", .)
    if (length(summitFiles) > 1) {
      nonOverlapES$Reproducibility <- rowSums(overlapMat)
      nonOverlapES$ReproducibilityPercent <- round(rowSums(overlapMat)/ncol(overlapMat), 
                                                   3)
      n <- length(summitFiles)
      minRep <- eval(parse(text = reproducibility))
      if (!is.numeric(minRep)) {
        stop("Error reproducibility not numeric when evaluated!")
      }
      idxPass <- which(nonOverlapES$Reproducibility >= 
                         minRep)
      nonOverlapPassES <- nonOverlapES[idxPass]
    }
    else {
      nonOverlapES$Reproducibility <- rep(NA, length(nonOverlapES))
      nonOverlapPassES <- nonOverlapES
    }
    nonOverlapPassES$groupScoreQuantile <- round(.getQuantiles(nonOverlapPassES$replicateScoreQuantile), 
                                                 3)
    mcols(nonOverlapPassES) <- mcols(nonOverlapPassES)[, 
                                                       c("score", "replicateScoreQuantile", "groupScoreQuantile", 
                                                         "Reproducibility", "GroupReplicate")]
    nonOverlapPassES
  })
  return(nonOverlapPassES)
}
test_Peak <- paste0('test_Peak/',list.files('test_Peak/'))
test_peak_name <- stringr::str_extract(list.files('test_Peak/'),'.*(?=_rmdup_summits)')
all_ReproduciblePeaks <- indentifyReproduciblePeaks(summitFiles = test_Peak,summitNames = test_peak_name)
all_ReproduciblePeaks_gtf <- generate_gtf(all_ReproduciblePeaks)
write.table(all_ReproduciblePeaks_gtf,file = 'all_ReproduciblePeaks.gtf',sep = "\t",col.names = F,row.names = F,quote = F)

#featureCounts -p -a all_ReproduciblePeaks.gtf -t exon -g gene_id -M -O -o all_peak_count.txt -T 20 ../rmdup_bam/