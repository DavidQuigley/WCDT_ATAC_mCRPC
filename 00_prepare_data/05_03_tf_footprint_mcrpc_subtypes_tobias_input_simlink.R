###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### LOAD LIBRARIES ---
library("stringr")
library("data.table")

### DEFINE PATH ---
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/analysis/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

dir.footprints <- file.path(dir.reproduce_data, "tf_footprints/mcrpc_subtypes")
dir.bam <- file.path(dir.reproduce_data, "tf_footprints/mcrpc_subtypes/bam")
dir.bed <- file.path(dir.reproduce_data, "tf_footprints/mcrpc_subtypes/bed")

dir.scripts <- file.path(dir.reproduce, "scripts_temp")  # REPLACE PATH
dir.logs <- tempdir() # REPLACE PATH

### DEFINE FILES ---
file.bed <- file.path(dir.bed, "mcrpc_subtypes_atacseq_peaks_uncorrected_merged.bed")

### GROUPS ---
subtypes <- c("ARpNEn","ARpNEp","ARlNEn","ARnNEp","ARnNEn")

### FUNCTION: getInputFiles() ---
getInputFiles <- function(dir.bam, dir.data, file.bed, subtype){
    # GET BAM FILES ----
    file.source.bam <- file.path(file.path(dir.bam, subtype), paste(subtype, ".bam", sep="") )
    file.source.bai <- file.path(file.path(dir.bam, subtype), paste(subtype, ".bam.bai", sep="") )
    file.target.bam <- file.path(dir.data, paste(subtype, ".bam", sep="") )
    file.target.bai <- file.path(dir.data, paste(subtype, ".bam.bai", sep="") )

    # GET BED FILE ---
    file.source.bed <- file.bed
    file.target.bed <- file.path(dir.data, "mcrpc_subtypes_atacseq_peaks_uncorrected_merged.bed")

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
for(subtype in subtypes){
    # OUTPUT DIR ---
    dir.subtype <- file.path(dir.footprints, subtype)
    dir.create(dir.subtype, showWarnings=FALSE)

    # GET INPUT FILES ---
    getInputFiles(dir.bam, dir.data=dir.subtype, file.bed, subtype)

    cat("PROCESSED:", subtype, "\n", sep="\t")
}


