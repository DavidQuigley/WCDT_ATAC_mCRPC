###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### LOAD LIBRARIES ---
library("stringr")
library("data.table")

### DEFINE PATH ---
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

dir.footprint <- file.path(dir.reproduce_data, "tf_footprints/mcrpc_samplewise")
dir.bed <- file.path(dir.footprints, "bed")

### DEFINE FILES ----
file.rds_atacseq_masterdata <- file.path(dir.reproduce_data, "wcdt_atacseq_masterdata.rds")
file.utility_functions <- file.path(dir.reproduce_scripts, "00_prepare_data/00_01_utility_functions.R")

### LOAD FUNCTIONS ---
source(file.utility_functions)

### CHROMOSOMES ---
chromosome <- paste("chr", c(1:22,"X","Y"), sep="")




### FUNCTION: mergeBEDregions() ---
mergeBEDregions <- function(list.peaks, sampleids){
    # GET NON-OVERLAPING REGIONS ---
    gr_nonoverlaping <- getNonOverlapPeaks(list.peaks=list.peaks, sampleids=sampleids)
    gr_nonoverlaping <- sort(gr_nonoverlaping)

    # GET DATA FRAME ---
    df.bed <- as.data.frame(gr_nonoverlaping)
    df.bed <- subset(df.bed, select=c("seqnames","start","end"))

    return(df.bed)
}



######################################################################
### LOAD ATAC-SEQ MASTER DATA ---
list.rds_atacseq_masterdata <- readRDS(file=file.rds_atacseq_masterdata)

### GET METADATA ---
metadata <- list.rds_atacseq_masterdata$metadata_mcrpc

### GET ATAC-SEQ PEAKS PER SAMPLE ---
list.peaks <- list.rds_atacseq_masterdata$atacseq_peaks_uncorrected



######################################################################
dat <- mergeBEDregions(list.peaks, sampleids=metadata$Sample_Name)

### WRITE OUTPUT ----
file.bed <- file.path(dir.bed, "mcrpc_samplewise_atacseq_peaks_uncorrected_merged.bed")
write.table(dat, file.bed, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)



######################################################################
### EXTRACT BED PER SAMPLE ----
for(i in 1:length(list.peaks)){
    sample_name <- names(list.peaks)[i]
    dir.output <- file.path(dir.footprints,  sample_name)
    file.bed <- file.path(dir.output, paste(sample_name, ".bed", sep=""))

    bed <- list.peaks[[sample_name]]
    bed$SampleID <- NULL
    bed$group <- "mcrpc"

    write.table(bed, file.bed, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

    cat("BED PROCESSED:", sample_name, "\n", sep="\t")
}


