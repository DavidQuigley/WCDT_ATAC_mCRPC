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

dir.bam <- file.path(dir.reproduce_data, "tf_footprints/mcrpc_subtypes/bam")
dir.bed <- file.path(dir.reproduce_data, "tf_footprints/mcrpc_subtypes/bed")

dir.scripts <- file.path(dir.reproduce, "scripts_temp")  # REPLACE PATH
dir.logs <- tempdir() # REPLACE PATH

### DEFINE FILES ----
file.rds_atacseq_masterdata <- file.path(dir.reproduce_data, "wcdt_atacseq_masterdata.rds")

### SUBTYPES ---
subtypes <- c("ARpNEn","ARpNEp","ARlNEn","ARnNEp","ARnNEn")
group_ids <- c("AR+_NE-","AR+_NE+","ARL_NE-","AR-_NE+","AR-_NE-")

### FUNCTION: get.SGEcmd()---
get.SGEcmd <- function(jobid, dir.logs){
    jobid <- jobid
    line01 <- "#!/bin/bash"

    ext_1 <- "MAX_MEM_TOTAL=30G"
    ext_2 <- "export SOFTWARE_BASE=/opt/software"
    ext_3 <- "export PICARD_PATH=${SOFTWARE_BASE}/picard/picard_2-7-1.jar"
    ext_4 <- "export PICARD=\"java -Xmx${MAX_MEM_TOTAL} -jar ${PICARD_PATH}\""

    line02 <- "#$ -pe smp 5"
    line03 <- "#$ -V"
    line04 <- "#$ -R y"
    line05 <- "#$ -l mem_free=30G"
    line06 <- paste("#$ -N", jobid, sep=" ")
    line07 <- paste("#$ -o", file.path(dir.logs, paste(jobid, ".out", sep="")), sep=" ")
    line08 <- paste("#$ -e", file.path(dir.logs, paste(jobid, ".err", sep="")), sep=" ")

    # COMPILE CMD ---
    cmd <- list(line01, ext_1, ext_2, ext_3, ext_4, line02, line03, line04, line05, line06, line07, line08)

    return(cmd)    
}


### FUNCTION: getMergeBAMs() ---
getMergeBAMs <- function(sampleids, dir.process, dir.output, subtype){
    # GET BAM FILES ---
    bamfile_suffix <- ".bowtie.sorted.nodup.bam"
    files.input_bam <- sapply(sampleids, function(x) {
                            paste(dir.process, "/", sprintf("%s/alignment/%s", x, x), bamfile_suffix, sep="")
                        })

    # OUTPUT BAM FILE ---
    file.output_bam <- file.path(dir.output, paste(subtype, ".bam", sep="") )

    # COMMANDS ---
    ln_1 <- paste("${PICARD} MergeSamFiles", "\\", sep=" ")
    ln_2 <- paste(paste("  I", files.input_bam, sep="="), "\\", sep=" ")
    ln_3 <- paste(paste("  O", file.output_bam, sep="="), "\\", sep=" ")
    ln_4 <- paste("  TMP_DIR=/scratch", "\\", sep=" ")
    ln_5 <- paste("  USE_THREADING=true", "\\", sep=" ")
    ln_6 <- "  CREATE_INDEX=true"

    # COMPILE TO LIST ---
    list.cmd <- list(ln_1, ln_2, ln_3, ln_4, ln_5, ln_6)

    return(list.cmd)
}

### FUNCTION: getMergeBAMsCMD() ---
getMergeBAMsCMD <- function(sampleids, dir.process, dir.logs, dir.output, subtype){
    # GET COMMANDS ---
    list.cmd_sge <- get.SGEcmd(jobid=paste("j", subtype, sep="_"), dir.logs)
    list.cmd_bam <- getMergeBAMs(sampleids, dir.process, dir.output, subtype)

    # MERGE COMMANDS ---
    cmd <- c(unlist(list.cmd_sge), unlist(list.cmd_bam) )

    return(cmd)
}



######################################################################
### LOAD ATAC-SEQ MASTER DATA ---
list.rds_atacseq_masterdata <- readRDS(file=file.rds_atacseq_masterdata)
metadata <- list.rds_atacseq_masterdata$metadata_mcrpc

### GET MERGE BAM COMMANDS ---
list.sh <- list()
for(i in 1:length(subtypes)){
    subtype <- subtypes[i]
    group_id <- group_ids[i]

    # GET COMMANDS ---
    sampleids <- metadata$Sample_Name[which(metadata$mCRPC_Subtype == group_id)]
    cmd <- getMergeBAMsCMD(sampleids, dir.process, dir.logs, dir.output=dir.bam, subtype=subtype)

    # WRITE QSUB SCRIPT --
    file.qsub <- file.path(dir.scripts, paste("qsub_mergebams", subtype, ".sh", sep=""))
    write.table(cmd, file.qsub, row.names=FALSE, col.names=FALSE, quote=FALSE)

    # STORE SCRIPT ---
    list.sh[[subtype]] <- paste("qsub", file.qsub, sep=" ")

    cat("PROCESSED:", subtype, "\n", sep=" ")
}


### WRITE MAIN RUN SCRIPT ---
file.sh <- file.path(dir.scripts, "main_mergebams.sh")
write.table(unlist(list.sh), file.sh, row.names=FALSE, col.names=FALSE, quote=FALSE)

cat("FILE GENERATED:", file.sh, "\n", sep=" ")
