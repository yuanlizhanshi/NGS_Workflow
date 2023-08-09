suppressWarnings(library(optparse))
option_list <- list(
  make_option(
    c("-s", "--species"),
    type="character", default=NULL,
    help = "Choose the species:[default]"
  )
)

args <- parse_args(OptionParser(option_list=option_list))

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

if (args$species == 'hs') {
  species <- "Homo sapiens"
}else if(args$species == 'mm'){
  species <- "Mus musculus"
}else{
  stop('Only support for hs and mm now')
}

fetch_trRNA_from_NCBI(species = species)
