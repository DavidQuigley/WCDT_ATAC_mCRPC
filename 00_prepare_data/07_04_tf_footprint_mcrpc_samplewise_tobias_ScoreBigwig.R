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

dir.scripts <- file.path(dir.reproduce, "scripts_temp")  # REPLACE PATH
dir.logs <- tempdir() # REPLACE PATH

### DEFINE FILES ----
file.rds_atacseq_masterdata <- file.path(dir.reproduce_data, "wcdt_atacseq_masterdata.rds")
file.bed <- file.path(dir.bed, "mcrpc_samplewise_atacseq_peaks_uncorrected_merged.bed")
file.ref_genome <- file.path("/opt/reference/Homo_sapiens/NCBI/GRCh38Decoy/Sequence/WholeGenomeFasta/genome.fa") # HUMAN REF GENOME HG38
file.blacklist_region <- file.path(dir.reproduce, "reference/GRCh38_unified_blacklist_2020_05_05.bed") # ENCODE BLACKLIST HG38
file.motif_db <- file.path(dir.reproduce, "reference/JASPAR2018_CORE_vertebrates_non-redundant.meme")
file.header_peaks <- file.path(dir.bed, "peaks_header_annotation.txt")

### CHROMOSOMES ---
chromosome <- paste("chr", c(1:22,"X","Y"), sep="")


######################################################################
### LOAD ATAC-SEQ MASTER DATA ---
list.rds_atacseq_masterdata <- readRDS(file=file.rds_atacseq_masterdata)

### GET METADATA ---
metadata <- list.rds_atacseq_masterdata$metadata_mcrpc
sample_names <- metadata$Sample_Name

### FUNCTION: get.SGEcmd()---
get.SGEcmd <- function(jobid, dir.logs){
    jobid <- jobid
    line01 <- "#!/bin/bash"
    line02 <- "#$ -pe smp 5"
    line03 <- "#$ -V"
    line04 <- "#$ -R y"
    line05 <- "#$ -l mem_free=10G"
    line06 <- paste("#$ -N", jobid, sep=" ")
    line07 <- paste("#$ -o", file.path(dir.logs, paste(jobid, ".out", sep="")), sep=" ")
    line08 <- paste("#$ -e", file.path(dir.logs, paste(jobid, ".err", sep="")), sep=" ")

    # COMPILE CMD ---
    cmd <- list(line01, line02, line03, line04, line05, line06, line07, line08)

    return(cmd)    
}


### FUNCTION: getScoreBigwig() ---
getScoreBigwig <- function(file.bigwig_corrected, file.bed_peaks, file.bigwig_footprints){
    # DEFINE VARS ---
    ln_1 <-  paste("FILE_BIGWIG_CORRECTED", file.bigwig_corrected, sep="=")
    ln_2 <-  paste("FILE_BIGWIG_FOOTPRINTS", file.bigwig_footprints, sep="=")
    ln_3 <-  paste("FILE_BED_PEAKS", file.bed_peaks, sep="=")

    # TOBIAS FUNCTION ---
    line1 <- paste("TOBIAS ScoreBigwig", "\\", sep=" ")
    line2 <- paste("  --signal ${FILE_BIGWIG_CORRECTED}", "\\", sep=" ")
    line3 <- paste("  --output ${FILE_BIGWIG_FOOTPRINTS}", "\\", sep=" ")
    line4 <- paste("  --regions ${FILE_BED_PEAKS}", "\\", sep=" ")
    line5 <- paste("  --score footprint", "\\", sep=" ")
    line6 <- paste("  --cores 40", "\\", sep=" ")
    line7 <- paste("  --verbosity 5", sep=" ")

    # COMPILE COMMAND ---
    list.cmd <- list(ln_1, ln_2, ln_3,
                    line1, line2, line3, line4, line5, line6, line7)

    return(list.cmd)
}

### GET FILES ---
files.bigwig_corrected <- sapply(sample_names, function(x) { paste(dir.footprint, "/", sprintf("%s/%s", x, x), "_corrected.bw", sep="") })
files.bigwig_footprints <- sapply(sample_names, function(x) { paste(dir.footprint, "/", sprintf("%s/%s", x, x), "_footprints.bw", sep="") })
files.bed <- sapply(sample_names, function(x) { paste(dir.footprint, "/", sprintf("%s/", x), "mcrpc_samplewise_atacseq_peaks_uncorrected_merged.bed", sep="") })


### GET COMMAND BY GROUPIDS ---
list.sh <- list()
for(sample_name in sample_names){
    # GET SOURCE FILE ---
    file.bigwig_corrected <- files.bigwig_corrected[sample_name]
    file.bigwig_footprints <- files.bigwig_footprints[sample_name]
    file.bed_peaks <- files.bed[sample_name]
    dir.output <- file.path(dir.footprint, sample_name)

    # GET COMMANDS ---
    list.cmd_sge <- get.SGEcmd(jobid=paste("ts", sample_name, sep="_"), dir.logs)
    list.cmd_tobias <- getScoreBigwig(file.bigwig_corrected, file.bed_peaks, file.bigwig_footprints)

    # MERGE COMMANDS ---
    cmd <- c(unlist(list.cmd_sge), unlist(list.cmd_tobias))

    # WRITE QSUB SCRIPT --
    file.qsub <- file.path(dir.scripts, paste("qsub_tobias_score_", sample_name, ".sh", sep=""))
    write.table(cmd, file.qsub, row.names=FALSE, col.names=FALSE, quote=FALSE)

    list.sh[[sample_name]] <- paste("qsub", file.qsub, sep=" ")

    cat("PROCESSED:", sample_name, "\n", sep=" ")
}

### WRITE MAIN RUN SCRIPT ---
file.sh <- file.path(dir.scripts, "main_tobias_score.sh")
write.table(unlist(list.sh), file.sh, row.names=FALSE, col.names=FALSE, quote=FALSE)

cat("FILE GENERATED:", file.sh, "\n", sep=" ")
