###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### DEFINE LIBRARIES ---
library("stringr")

### DEFINE PATH ---
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/analysis/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

dir.tfpeaks2gene <- file.path(dir.reproduce_data, "tfpeaks2gene")
dir.footprints <- file.path(dir.reproduce_data, "tf_footprints/mcrpc_subtypes/data_footprints") 
dir.tf2genes_intermediate <- file.path(dir.reproduce_data, "tf_footprints/mcrpc_subtypes/tf2genes_intermediate") #CREATE

dir.scripts <- file.path(dir.reproduce, "scripts_temp")  # REPLACE PATH
dir.logs <- tempdir() # REPLACE PATH


### DEFINE FILE ---
file.peaks2gene <- file.path(dir.tfpeaks2gene, "wcdt_tfpeaks2genes_correlated_r_0p4_pval_0p05_distTSS_500kb.tsv.gz")
file.rds_tf2peaks <- file.path(dir.tfpeaks2gene, "tf2peaks_mcrpc_subtypes.rds")
file.script <- file.path(dir.reproduce_scripts, "11_01_tf2genes_links_function.R")

### FUNCTION: get.SGEcmd()---
get.SGEcmd <- function(jobid, dir.logs){
    jobid <- jobid
    line01 <- "#!/bin/bash"
    line02 <- "#$ -pe smp 1"
    line03 <- "#$ -V"
    line04 <- "#$ -R y"
    line05 <- "#$ -l mem_free=1G"
    line06 <- paste("#$ -N", jobid, sep=" ")
    line07 <- paste("#$ -o", file.path(dir.logs, paste(jobid, ".out", sep="")), sep=" ")
    line08 <- paste("#$ -e", file.path(dir.logs, paste(jobid, ".err", sep="")), sep=" ")

    # COMPILE CMD ---
    cmd <- list(line01, line02, line03, line04, line05, line06, line07, line08)

    return(cmd)    
}



### LOAD TF2PEAKS DATA ---
list.tf2peaks <- readRDS(file=file.rds_tf2peaks)
motif_ids <- names(list.tf2peaks)


### GET TF2GENES FOR EACH MOTIF ---
list.sh <- list()
for(k in 1:length(motif_ids)){
    motif_prefix <- motif_ids[k]
    file.tf2peaks <- file.rds_tf2peaks
    file.output <- file.path(dir.tf2genes_intermediate, paste("tf2genes_", k, ".tsv", sep=""))

    # GET COMMANDS ---
    list.cmd_sge <- get.SGEcmd(jobid=paste("t2g", k, sep="_"), dir.logs)
    cmd_r <- paste("Rscript", file.script, 
                        "-t", file.tf2peaks, 
                        "-p", file.peaks2gene,
                        "-m", motif_prefix, 
                        "-o", file.output, 
                        sep=" ")

    # MERGE COMMANDS ---
    cmd <- c(unlist(list.cmd_sge), cmd_r)

    # WRITE QSUB SCRIPT --
    file.qsub <- file.path(dir.scripts, paste("qsub_", k, ".sh", sep=""))
    write.table(cmd, file.qsub, row.names=FALSE, col.names=FALSE, quote=FALSE)

    list.sh[[k]] <- paste("qsub", file.qsub, sep=" ")

    cat("PROCESSED:", k, "\n", sep=" ")
}

### WRITE MAIN RUN SCRIPT ---
file.sh <- file.path(dir.scripts, "main.sh")
write.table(unlist(list.sh), file.sh, row.names=FALSE, col.names=FALSE, quote=FALSE)

cat("FILE GENERATED:", file.sh, "\n", sep=" ")
