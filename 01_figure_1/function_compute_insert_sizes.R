###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### LOAD LIBRARIES ---
suppressPackageStartupMessages(require("stringr", quietly=TRUE, warn.conflicts=FALSE))
suppressPackageStartupMessages(require("dplyr", quietly=TRUE, warn.conflicts=FALSE))
suppressPackageStartupMessages(require("GenomicRanges", quietly=TRUE, warn.conflicts=FALSE))
suppressPackageStartupMessages(require("GenomicAlignments", quietly=TRUE, warn.conflicts=FALSE))
suppressPackageStartupMessages(require("Rsamtools", quietly=TRUE, warn.conflicts=FALSE))
suppressPackageStartupMessages(require("optparse", quietly=TRUE, warn.conflicts=FALSE))

### PASS ARGUMENTS ---
list.option <- list(
  optparse::make_option(c("-b", "--file_bam"), type = "character", help = "Path to BAM file [Required]"),
  optparse::make_option(c("-o", "--file_output"), type = "character", help = "Path to output file [Required]"),
  optparse::make_option(c("-s", "--sampleid"), type = "character", help = "SampleID [Required]")
)

parseobj <- optparse::OptionParser(option_list=list.option, usage="usage: Rscript %prog[options]")
opt <- parse_args(parseobj)

### PASS OBJECT VALUES ---
file.bam <- opt$file_bam
file.output <- opt$file_output
sampleid <- opt$sampleid

### FUNCTION: getInsertSizesbySampleID() ---
getInsertSizesbySampleID <- function(file.bam, sampleid){
    chromosomes <- paste("chr", c(1:22,"X","Y"), sep="")

    # SET BAM SCANNING PARAMETERS ---
    param <- Rsamtools::ScanBamParam(mapqFilter = 1, flag = scanBamFlag(isPaired = TRUE, isProperPair = TRUE), what = c("qname", "mapq", "isize"))

    # READ BAM FILES ---
    g_align <- GenomicAlignments::readGAlignmentPairs(file.bam, param=param)

    # GET FIRST PAIR-END READS ---
    g_align_reads1 <- GenomicAlignments::first(g_align)

    # ONLY KEEP READS IN chr1-chrY ----
    g_align_reads1_filtered <- g_align_reads1[GenomicAlignments::seqnames(seqinfo(g_align_reads1)) %in% chromosomes]

    # GET INSERT SIZE ---
    insertSizes <- abs(GenomicRanges::elementMetadata(g_align_reads1_filtered)$isize)

    # GET INSERT SIZE FREQUENCY TABLE ---
    dm <- data.frame(insertSizes) %>%
            dplyr::count(insertSizes, sort=FALSE, name="Count") %>% 
            dplyr::mutate(CountLog= log10(Count)) %>% 
            dplyr::mutate(Proportion=(Count/sum(Count))*100)
    colnames(dm) <- c("InsertSizes","Count","CountLog","Proportion")

    ## ADD SAMPLEID ---
    dm$SampleID <- sampleid

    return(dm)
}

### GET INSERT SIZES ---
dm <- getInsertSizesbySampleID(file.bam, sampleid)

### WRITE OUPTUT ---
write.table(dm, file.output, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
