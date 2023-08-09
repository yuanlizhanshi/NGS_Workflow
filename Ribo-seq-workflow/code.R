library(data.table)
library(tidyverse)

##Doownload srr
df <- fread('sra_explorer_sra_download.sh')
srr <- str_extract(df$V5,'SRR\\d+')
srr2 <- paste0('https://sra-pub-run-odp.s3.amazonaws.com/sra/',srr,'/',srr)
write_lines(srr2,file = 'srr2.txt')
##Download rRNA and tRNA
fetch_trRNA_from_NCBI <- function(species = "Homo sapiens",
                                  output_file = NULL){
  # search
  rs <- rentrez::entrez_search(db = "nucleotide",
                               term = paste("(",species,"[ORGN] AND rRNA[FILT]) OR (",
                                            species,"[ORGN] AND tRNA[FILT])",sep = ""),
                               retmax = 20000,
                               use_history = TRUE)

  # download sequence
  if(is.null(output_file)){
    output_file <- paste(gsub(" ",replacement = "_",species),"_trRNA.fa",sep = "")
    file.create(output_file)
  }else{
    file.create(output_file)
  }

  ids = rs$ids
  fa <- rentrez::entrez_fetch(db = "nucleotide",rettype = "fasta",
                              web_history = rs$web_history)
  write(fa,file = output_file,append = TRUE)

  raw_fa <- Biostrings::readDNAStringSet(output_file)
  Biostrings::writeXStringSet(raw_fa,filepath = output_file)
  cat("Download tRNA and rRNA sequences from NCBI has finished!")
}
fetch_trRNA_from_NCBI(species = 'Mus musculus')
setClass(Class = "GenomeGTF",
         slots = c(gtf = "GRanges",
                   gtfPath = "character",
                   genome = "DNAStringSet",
                   genomePath = "character",
                   representTrans = "data.frame",
                   intron = "data.frame",
                   ranges = "GroupedIRanges"),
         prototype = list(gtf = NULL,
                          genome = NULL,
                          gtfPath = NULL,
                          genomePath = NULL),
         contains = c("DNAStringSet"))

loadGenomeGTF <- function(gtfPath = NULL,
                          genomePath = NULL,
                          format = "gtf",
                          filterProtein = FALSE){
  
  
  setClass(Class = "GenomeGTF",
           slots = c(gtf = "GRanges",
                     gtfPath = "character",
                     genome = "DNAStringSet",
                     genomePath = "character",
                     representTrans = "data.frame",
                     intron = "data.frame",
                     ranges = "GroupedIRanges"),
           prototype = list(gtf = NULL,
                            genome = NULL,
                            gtfPath = NULL,
                            genomePath = NULL),
           contains = c("DNAStringSet"))
  # load gtf
  gtf_input <- rtracklayer::import.gff(gtfPath,format = format)
  # as.data.frame()
  
  # load genome sequences
  if(!is.null(genomePath)){
    myFASTA <- Biostrings::readDNAStringSet(genomePath,format = "fasta")
  }else{
    genomePath = ""
    myFASTA = Biostrings::DNAStringSet(NULL)
  }
  
  # whether get protein
  if(filterProtein == TRUE){
    if("ccds_id" %in% colnames(data.frame(gtf_input))){
      protein <- gtf_input[which(gtf_input$ccds_id != "NA"),]
      gtf <- gtf_input[which(gtf_input$transcript_id %in% unique(protein$transcript_id)),]
    }else{
      protein <- gtf_input[which(gtf_input$type == "CDS"),]
      gtf <- gtf_input[which(gtf_input$gene_id %in% unique(protein$gene_id)),]
    }
    
  }else{
    gtf <- gtf_input
  }

  object <-
    methods::new("GenomeGTF",
                 gtf = gtf,
                 gtfPath = gtfPath,
                 genome =  myFASTA,
                 genomePath = genomePath)
  
  return(object)
}

filterID <- function(object,
                     geneName = NULL,geneId = NULL,transId = NULL){
  # load gtf
  gtf <- as.data.frame(object@gtf)
  
  if(!is.null(geneName) & is.null(geneId) & is.null(transId)){
    ginfo <- gtf[which(gtf$gene_name %in% geneName),]
  }else if(is.null(geneName) & !is.null(geneId) & is.null(transId)){
    ginfo <- gtf[which(gtf$gene_id %in% geneId),]
  }else if(is.null(geneName) & is.null(geneId) & !is.null(transId)){
    ginfo <- gtf[which(gtf$transcript_id %in% transId),]
  }else if(!is.null(geneName) & is.null(geneId) & !is.null(transId)){
    ginfo <- gtf[which(gtf$gene_name == geneName & gtf$transcript_id == transId),]
  }else if(is.null(geneName) & !is.null(geneId) & !is.null(transId)){
    ginfo <- gtf[which(gtf$gene_id == geneId & gtf$transcript_id == transId),]
  }else{
    message("Please choose again.")
  }
  
  # output
  return(ginfo)
}
getTransInfo <- function(object,
                   geneName = NULL,geneId = NULL,transId = NULL,
                   selecType = c("lcds","lt"),
                   topN = 1,sep = "|",
                   filterGene = FALSE){

            # load GTF
            ginfo <- filterID(object = object,geneName = geneName,geneId = geneId,transId = transId)
            
            # whether filter gene not in genome file chromosomes
            if(filterGene == TRUE){
              chrName <- names(object@genome)
              ginfo <- ginfo %>%
                dplyr::filter(seqnames %in% chrName)
            }
            # ===============================================================================
            # recode
            # prepare test type
            test <- ginfo[1:4,] %>% dplyr::mutate(gene_id = "test",
                                                  transcript_id = "test")
            test$type <- c("exon","5UTR","CDS","3UTR")
            ginfo <- rbind(test,ginfo)
            
            # get info
            leninfo <- ginfo[which(ginfo$type %in% c("exon","5UTR","five_prime_utr","CDS","3UTR","three_prime_utr")),] %>%
              dplyr::mutate(type = dplyr::if_else(type == "five_prime_utr","5UTR",type)) %>%
              dplyr::mutate(type = dplyr::if_else(type == "three_prime_utr","3UTR",type)) %>%
              dplyr::group_by(gene,gene_id,transcript_id,type) %>%
              dplyr::summarise(typelen = sum(width)) %>%
              tidyr::spread(.,type,typelen,fill = 0) %>%
              dplyr::mutate(gtype = dplyr::if_else(`5UTR` > 0 | `CDS` > 0,"CD","NC")) %>%
              dplyr::mutate(cdsst = dplyr::if_else(`gtype` == "CD",`5UTR` + 1,1),
                            cdsed = dplyr::if_else(`gtype` == "CD",`cdsst` + `CDS`,exon)) %>%
              dplyr::mutate(cdsed = dplyr::if_else(cdsed > exon,cdsed - 3,cdsed)) %>%
              dplyr::mutate(tname = paste(gene,gene_id,transcript_id,cdsst,cdsed,exon,gtype,sep = sep)) %>%
              dplyr::filter(gene_id != "test")
            
            # ===============================================================================
            # filter
            if(topN == 0){
              final.res <- leninfo
            }else{
              # choose type
              if(selecType == "lt"){
                final.res <- leninfo %>%
                  dplyr::group_by(gene,gene_id) %>%
                  dplyr::arrange(dplyr::desc(exon),dplyr::desc(CDS)) %>%
                  dplyr::slice_head(n = as.numeric(topN))
              }else if(selecType == "lcds"){
                final.res <- leninfo %>%
                  dplyr::group_by(gene,gene_id) %>%
                  dplyr::arrange(dplyr::desc(CDS),dplyr::desc(exon)) %>%
                  dplyr::slice_head(n = as.numeric(topN))
              }else{
                message("Please choose 'lt' or 'lcds'!")
              }
            }
            
            return(final.res)
          }

gtf <- loadGenomeGTF(gtfPath = "index/mm10.gtf",filterProtein = T)
gene <- unique(gtf@gtf$gene_id)
rt <- getTransInfo(object = gtf,geneId = gene,
                   selecType = "lcds",topN = 1)
rt <- rt %>% 
  mutate(seq_id = paste(gene,gene_id,transcript_id,
                        `5UTR` + 1,`5UTR` + CDS,exon,sep = "|"))
fwrite(rt,file = 'seq_info.txt',sep = '\t')


#Rscript.exe Download_trRNA.R -s mm
#bowtie2_build Mus_musculus_trRNA.fa mm_trRNA
#Rscript.exe get_seqinfo.R -g index/mm10.gtf
#python getSeq.py index/mm10.gtf index/mm10.fa seq_info.txt
#hisat2-build longest_transcript.fasta mm_longest_transcript



library(RiboProfiler)
pre_longest_trans_info(gtf_file = "index/mm10.gtf",
                       out_file = "longest_info.txt")



export_gtf <- function(gtf_file = NULL){
  gtf <- rtracklayer::import.gff(gtf_file, format = "gtf")
  gtf <- data.frame(gtf)
  colnames(gtf)[14] <- 'gene_name'
  gtf_plus <- dplyr::arrange(dplyr::filter(gtf, 
                                           strand %in% "+"), seqnames, gene_name, gene_id, 
                             transcript_id, type, start, end)
  gtf_neg <- dplyr::arrange(dplyr::filter(gtf, 
                                          strand %in% "-"), seqnames, gene_name, gene_id, 
                            transcript_id, type, dplyr::desc(start), dplyr::desc(end))
  sorted_gtf <- rbind(gtf_plus, gtf_neg)
  output_name = paste(gtf_file, ".sorted.gtf", sep = "")
  rtracklayer::export.gff(sorted_gtf, con = output_name, 
                          format = "gtf")
}
export_gtf('index/mm10.gtf')
