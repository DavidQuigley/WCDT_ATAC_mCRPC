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

dir.input <- file.path(dir.reproduce_data, "tfpeaks2gene/input")
dir.corr <- file.path(dir.reproduce_data, "tfpeaks2gene/correlation") 

dir.scripts <- file.path(dir.reproduce, "scripts_temp")  # REPLACE PATH
dir.logs <- tempdir() # REPLACE PATH

### DEFINE FILES ---
file.script <- file.path(dir.reproduce_scripts, "08_03_tfpeaks2genes_correlation_function.R")


### CHROMOSOME ---
chromosomes <- paste("chr", c(1:22, "X","Y"), sep="")

### FUNCTION: get.SGEcmd()---
get.SGEcmd <- function(jobid, dir.logs){
    jobid <- jobid
    line01 <- "#!/bin/bash"
    line02 <- "#$ -pe smp 2"
    line03 <- "#$ -V"
    line04 <- "#$ -R y"
    line05 <- "#$ -l mem_free=6G"
    line06 <- paste("#$ -N", jobid, sep=" ")
    line07 <- paste("#$ -o", file.path(dir.logs, paste(jobid, ".out", sep="")), sep=" ")
    line08 <- paste("#$ -e", file.path(dir.logs, paste(jobid, ".err", sep="")), sep=" ")

    # COMPILE CMD ---
    cmd <- list(line01, line02, line03, line04, line05, line06, line07, line08)

    return(cmd)    
}



### GET PEAKS2GENE CORRELATION BY CHROMOSOME ---
list.sh <- list()
for(chr in chromosomes){
    file.atacseq <- file.path(dir.input, paste("atacseq_peaks2gene_input_", chr,  ".tsv.gz", sep=""))
    file.rnaseq <- file.path(dir.input, paste("rnaseq_peaks2gene_input_", chr,  ".tsv.gz", sep=""))
    file.output <- file.path(dir.corr, paste("peaks2gene_corr_", chr,  ".tsv", sep=""))

    # GET COMMANDS ---
    list.cmd_sge <- get.SGEcmd(jobid=paste("p2g", chr, sep="_"), dir.logs)
    cmd_r <- paste("Rscript", file.script, 
                        "-a", file.atacseq, 
                        "-r", file.rnaseq, 
                        "-o", file.output, sep=" ")

    # MERGE COMMANDS ---
    cmd <- c(unlist(list.cmd_sge), cmd_r)

    # WRITE QSUB SCRIPT --
    file.qsub <- file.path(dir.scripts, paste("qsub_p2g_", chr, ".sh", sep=""))
    write.table(cmd, file.qsub, row.names=FALSE, col.names=FALSE, quote=FALSE)

    list.sh[[chr]] <- paste("qsub", file.qsub, sep=" ")

    cat("PROCESSED:", chr, "\n", sep=" ")
}

### WRITE MAIN RUN SCRIPT ---
file.sh <- file.path(dir.scripts, "main_p2g.sh")
write.table(unlist(list.sh), file.sh, row.names=FALSE, col.names=FALSE, quote=FALSE)

cat("FILE GENERATED:", file.sh, "\n", sep=" ")

