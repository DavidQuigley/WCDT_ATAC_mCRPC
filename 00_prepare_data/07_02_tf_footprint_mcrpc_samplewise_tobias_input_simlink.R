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

dir.processed_atacseq <- file.path("/data1/datasets_1/human_prostate_WCDT/atac_solid/processed/2020_09_22/processed") # REPLACE THIS
dir.footprint <- file.path(dir.reproduce_data, "tf_footprints/mcrpc_samplewise")
dir.bed <- file.path(dir.footprints, "bed")

### DEFINE FILES ----
file.rds_atacseq_masterdata <- file.path(dir.reproduce_data, "wcdt_atacseq_masterdata.rds")
file.bed <- file.path(dir.bed, "mcrpc_samplewise_atacseq_peaks_uncorrected_merged.bed")

### CHROMOSOMES ---
chromosome <- paste("chr", c(1:22,"X","Y"), sep="")


######################################################################
### LOAD ATAC-SEQ MASTER DATA ---
list.rds_atacseq_masterdata <- readRDS(file=file.rds_atacseq_masterdata)

### GET METADATA ---
metadata <- list.rds_atacseq_masterdata$metadata_mcrpc

### GET BAM FILES ---
files.bam <- sapply(metadata$Sample_Name, function(x) { paste(dir.processed_atacseq, "/", sprintf("%s/alignment/%s", x, x), ".bowtie.sorted.nodup.bam", sep="") })





### FUNCTION: getInputFiles() ---
getInputFiles <- function(dir.data, files.bam, file.bed, sampleid){
    # GET BAM FILES ----
    file.source.bam <- files.bam[sampleid]
    file.source.bai <- paste(files.bam[sampleid], ".bai", sep="")
    file.target.bam <- file.path(dir.data, paste(sampleid, ".bam", sep="") )
    file.target.bai <- file.path(dir.data, paste(sampleid, ".bam.bai", sep="") )

    # GET BED FILE ---
    file.source.bed <- file.bed
    file.target.bed <- file.path(dir.data, "mcrpc_samplewise_atacseq_peaks_uncorrected_merged.bed")

    # CREATE SOFT LINK ---
    cmd1 <- paste("ln -s", file.source.bam, file.target.bam, sep=" ")
    cmd2 <- paste("ln -s", file.source.bai, file.target.bai, sep=" ")
    cmd3 <- paste("ln -s", file.source.bed, file.target.bed, sep=" ")
    
    # EXECUTE COMMAND ---
    system(cmd1)
    system(cmd2)    
    system(cmd3)    
}



### PREPARE SIMLINK ---
for(sample_name in metadata$Sample_Name){
    # OUTPUT DIR ---
    dir.sampleid <- file.path(dir.footprints, sample_name)
    dir.create(dir.sampleid, showWarnings=FALSE)

    # GET INPUT FILES ---
    getInputFiles(dir.data=dir.sampleid, files.bam, file.bed, sampleid=sample_name)

    cat("PROCESSED:", sample_name, "\n", sep="\t")
}

