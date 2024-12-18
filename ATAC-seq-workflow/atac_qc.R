suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ATACseqQC))
options(warn = -1)
option_list <- list(
    make_option(
        c("-b", "--bam"),
        type="character", default=NULL,
        help = "Input bam file"
    )
)
args <- parse_args(OptionParser(option_list=option_list))
bam_file =  args$bam
bamFiles.labels <- tools::file_path_sans_ext(basename(bam_file))
if (!dir.exists('fragment_distuibution')) {
    dir.create('fragment_distuibution')
}
out_pdf_name <- paste0('./fragment_distuibution/',bamFiles.labels,'.pdf')
pdf(out_pdf_name,width = 8,height = 6)
fragSize <- fragSizeDist(bamFile = bam_file,bamFiles.labels =  bamFiles.labels)
dev.off()
##equal to
##samtools view ATAC_f2q30_sorted.bam | awk '$9>0' | cut -f 9 | sort | uniq -c | sort -b -k2,2n | sed -e 's/^[ \t]*//' > fragment_length_count.txt
