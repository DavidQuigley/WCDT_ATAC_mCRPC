###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### DEFINE LIBRARIES ---
library("stringr")
library("Rsamtools")
library("GenomicAlignments")
library("GenomicRanges")
library("IRanges")

### DEFINE PATH ---
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

dir.processed_atacseq <- file.path("/data1/datasets_1/human_prostate_WCDT/atac_solid/processed/2020_09_22/processed") # REPLACE THIS BAM FILE PATH AS REQUIRED

### DEFINE FILES ----
file.rds_atacseq_masterdata <- file.path(dir.reproduce_data, "wcdt_atacseq_masterdata.rds")

### CHROMOSOMES ---
chromosomes <- paste("chr", c(1:22,"X","Y"), sep="")


####################################################################################################
### FUNCTION: compute_FRiP() ---
compute_FRiP <- function(file.bam, gr_peak, sampleid){
    # SET BAM SCANNING PARAMETERS ---
    param <- Rsamtools::ScanBamParam(flag = scanBamFlag(isPaired = TRUE, isProperPair = TRUE), what = c("qname", "pos", "qwidth", "mapq", "mpos", "isize"))

    # READ BAM FILES ---
    g_align <- GenomicAlignments::readGAlignmentPairs(file.bam, param=param)
    g_reads <- GenomicRanges::GRanges(g_align)

    # GET TOTAL NUMBER OF READS ---
    mapped_paired_reads <- length(GenomicRanges::ranges(g_reads))    

    # GET PEAK OVERLAPS ---
    overlap_counts <- IRanges::countOverlaps(g_reads, gr_peak)
    reads_in_peaks <- length(overlap_counts[which(overlap_counts > 0)])

    # COMPUTE FRiP RATIO ---
    frip <- reads_in_peaks/mapped_paired_reads

    # COMPILE OUTPUT ---
    df <- data.frame(SampleID=sampleid, 
                        Mapped_paired_reads=mapped_paired_reads, 
                        Reads_in_peaks=reads_in_peaks,
                        FRiP=frip)
    
    return(df)
}







####################################################################################################
### LOAD ATAC-SEQ MASTER DATA ---
list.rds_atacseq_masterdata <- readRDS(file=file.rds_atacseq_masterdata)

### GET METADATA ---
metadata <- list.rds_atacseq_masterdata$metadata_mcrpc
sampleids <- metadata$Sample_ID

### GET ATAC-SEQ PEAKS ---
list.peaks <- list.rds_atacseq_masterdata$atacseq_peaks
names(list.peaks) <- metadata$Sample_ID
list.gr_peaks <- lapply(list.peaks, function(x) GenomicRanges::makeGRangesFromDataFrame(x) )

### GET BAM FILES ---
files.bam <- sapply(metadata$Sample_Name, function(x) { paste(dir.processed_atacseq, "/", sprintf("%s/alignment/%s", x, x), ".bowtie.sorted.nodup.bam", sep="") })
names(files.bam) <- metadata$Sample_ID




####################################################################################################
##### GET FRiP SCORES FOR EACH SAMPLEID ---
list.dat <- list()
for(sampleid in sampleids){
    file.bam <- files.bam[sampleid]
    gr_peak <- list.gr_peaks[[sampleid]]
    
    # GET FRiP SCORES ---
    list.dat[[sampleid]] <- compute_FRiP(file.bam, gr_peak, sampleid)
}

### COMPILE OUTPUT ---
dat <- do.call(rbind.data.frame, list.dat)
rownames(dat) <- NULL





####################################################################################################
### WRITE OUTPUT ---
file.output <- ""
write.table(dat, file.output, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
