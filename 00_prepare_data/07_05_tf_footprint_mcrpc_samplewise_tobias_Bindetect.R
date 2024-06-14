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


### FUNCTION: getBindetect() ---
getBindetect <- function(dir.data, file.motif_db, file.ref_genome, file.bed_peaks, file.header_peaks, sampleid, group_name){
    # DEFINE FILES ---
    dir.sampleid <- file.path(dir.data, sampleid)
    dir.output_group <- file.path(dir.sampleid, "bindetect")
    file.bw <- file.path(dir.sampleid, paste(sampleid, "footprints.bw", sep="_")) 

    # DEFINE VARS ---
    ln_1 <-  paste("DIR_OUTPUT", dir.output_group, sep="=")
    ln_2 <-  paste("FILE_BIGWIG_1", file.bw, sep="=")
    ln_3 <-  paste("FILE_BED_PEAKS", file.bed_peaks, sep="=")
    ln_4 <-  paste("FILE_MOTIF_DB", file.motif_db, sep="=")
    ln_5 <-  paste("FILE_REF_GENOME", file.ref_genome, sep="=")
    ln_6 <-  paste("CONDITION_1", group_name, sep="=")
    ln_7 <-  paste("FILE_HEADER_PEAKS", file.header_peaks, sep="=")
    ln_8 <- "mkdir ${DIR_OUTPUT}"

    # TOBIAS FUNCTION ---
    line1 <- paste("TOBIAS BINDetect", "\\", sep=" ")
    line2 <- paste("  --signals ${FILE_BIGWIG_1}", "\\", sep=" ")
    line3 <- paste("  --peaks ${FILE_BED_PEAKS}", "\\", sep=" ")
    line4 <- paste("  --peak-header ${FILE_HEADER_PEAKS}", "\\", sep=" ")
    line5 <- paste("  --motifs ${FILE_MOTIF_DB}", "\\", sep=" ")
    line6 <- paste("  --genome ${FILE_REF_GENOME}", "\\", sep=" ")
    line7 <- paste("  --cond-names ${CONDITION_1}", "\\", sep=" ")
    line8 <- paste("  --naming \"id_name\" ", "\\", sep=" ")
    line9 <- paste("  --motif-pvalue 1e-4", "\\", sep=" ")
    line10 <- paste("  --bound-pvalue 0.001", "\\", sep=" ")
    line11 <- paste("  --outdir ${DIR_OUTPUT}", "\\", sep=" ")
    line12 <- paste("  --cores 50", "\\", sep=" ")
    line14 <- paste("  --verbosity 5", sep=" ")

    # COMPILE COMMAND ---
    list.cmd <- list(ln_1, ln_2, ln_3, ln_4, ln_5, ln_6, ln_7, ln_8,
                    line1, line2, line3, line4, line5, line6, line7, line8, line9, line10, line11, line12, line14)

    return(list.cmd)
}

### GET FILES ---
files.bed <- sapply(sample_names, function(x) { paste(dir.footprint, "/", sprintf("%s/%s", x,x), ".bed", sep="") })

### GET COMMAND BY GROUPIDS ---
list.sh <- list()
for(sample_name in sample_names){
    # GET SOURCE FILE ---
    file.bed_peaks <- files.bed[sample_name]

    # GET COMMANDS ---
    list.cmd_sge <- get.SGEcmd(jobid=paste("bd", sample_name, sep="_"), dir.logs)
    list.cmd_tobias <- getBindetect(dir.data=dir.footprint, file.motif_db, file.ref_genome, file.bed_peaks, file.header_peaks, sampleid=sample_name, group_name="mcrpc")

    # MERGE COMMANDS ---
    cmd <- c(unlist(list.cmd_sge), unlist(list.cmd_tobias))

    # WRITE QSUB SCRIPT --
    file.qsub <- file.path(dir.scripts, paste("qsub_tobias_bindetect_", sample_name, ".sh", sep=""))
    write.table(cmd, file.qsub, row.names=FALSE, col.names=FALSE, quote=FALSE)

    list.sh[[sample_name]] <- paste("qsub", file.qsub, sep=" ")

    cat("PROCESSED:", sample_name, "\n", sep=" ")
}

### WRITE MAIN RUN SCRIPT ---
file.sh <- file.path(dir.scripts, "main_tobias_bindetect.sh")
write.table(unlist(list.sh), file.sh, row.names=FALSE, col.names=FALSE, quote=FALSE)

cat("FILE GENERATED:", file.sh, "\n", sep=" ")
