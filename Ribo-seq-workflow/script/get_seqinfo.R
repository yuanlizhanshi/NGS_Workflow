suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Biostrings))
options(warn = -1)
option_list <- list(
  make_option(
    c("-g", "--gtf"),
    type="character", default=NULL,
    help = "Choose the gtf file"
  )
)

args <- parse_args(OptionParser(option_list=option_list))
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


gtf <- loadGenomeGTF(gtfPath = args$gtf,filterProtein = T)
gene <- unique(gtf@gtf$gene_id)
rt <- getTransInfo(object = gtf,geneId = gene,
                   selecType = "lcds",topN = 1)
rt <- rt %>%
  mutate(seq_id = paste(gene,gene_id,transcript_id,
                        `5UTR` + 1,`5UTR` + CDS,exon,sep = "|"))
fwrite(rt,file = 'seq_info.txt',sep = '\t')