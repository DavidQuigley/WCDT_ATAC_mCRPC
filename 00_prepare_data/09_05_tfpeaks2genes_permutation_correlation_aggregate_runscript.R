###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### DEFINE LIBRARIES ---
library("stringr")

### DEFINE PATH ---
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

dir.input <- file.path(dir.reproduce_data, "tfpeaks2gene/input")
dir.corr <- file.path(dir.reproduce_data, "tfpeaks2gene/correlation") 
dir.permute <- file.path(dir.reproduce_data, "tfpeaks2gene/permutation")
dir.permute_corr <- file.path(dir.reproduce_data, "tfpeaks2gene/permutation/perm_correlation")
dir.permute_agg <- file.path(dir.reproduce_data, "tfpeaks2gene/permutation/perm_correlation_aggregate")

dir.scripts <- file.path(dir.reproduce, "scripts_temp")  # REPLACE PATH
dir.logs <- tempdir() # REPLACE PATH

### DEFINE FILES ---
file.script <- file.path(dir.reproduce_scripts, "09_04_tfpeaks2genes_permutation_correlation_aggregate_function.R")
file.perm_rds <- file.path(dir.permute, "wcdt_rnaseq_permutation.rds")

### CHROMOSOME ---
chromosomes <- paste("chr", c(1:22, "X","Y"), sep="")



### FUNCTION: get.SGEcmd()---
get.SGEcmd <- function(jobid, dir.logs){
    jobid <- jobid
    line01 <- "#!/bin/bash"
    line02 <- "#$ -pe smp 10"
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



### GET PEAKS2GENE CORRELATION BY CHROMOSOME ---
list.sh <- list()
n <- 1000

for(chr in chromosomes)
    ### ---- CHECK FILE STATUS ----
    list.err <- list()
    for(k in 1:1000){
        key <- paste(chr, k, sep="_")
        file.dat <- file.path(dir.permute_corr, paste("perm_", key, ".tsv.gz", sep=""))
        
        if(!file.exists(file.dat)){
            list.err[[key]] <- file.dat
            cat("ERROR! FILE NOT FOUND:", file.dat, "\n", sep=" ")
        }
    }

    if(length(list.err) > 0 ){
        stop("SCRIPT TERMINATED PREMATURELY BECAUSE OF ERROR ABOVE!")
    }
    ### ---- CHECK FILE STATUS ----

    # DEFINE OUTPUT FILE ---
    file.output <- file.path(dir.permute_agg, paste("perm_agg", chr, ".tsv", sep=""))

    # GET COMMAND ---
    list.cmd_sge <- get.SGEcmd(jobid=paste("agg", chr, sep="_"), dir.logs)
    cmd_r <- paste("Rscript", file.script, 
                        "-d", dir.input, 
                        "-o", file.output, 
                        "-c", chr, sep=" ")

    # MERGE COMMANDS ---
    cmd <- c(unlist(list.cmd_sge), cmd_r)

    # WRITE QSUB SCRIPT --
    file.qsub <- file.path(dir.scripts, paste("qsub_", chr, ".sh", sep=""))
    write.table(cmd, file.qsub, row.names=FALSE, col.names=FALSE, quote=FALSE)

    list.sh[[chr]] <- paste("qsub", file.qsub, sep=" ")

    cat("PROCESSED:", chr, "\n", sep=" ")
}

### WRITE MAIN RUN SCRIPT ---
file.sh <- file.path(dir.scripts, "main.sh")
write.table(unlist(list.sh), file.sh, row.names=FALSE, col.names=FALSE, quote=FALSE)

cat("FILE GENERATED:", file.sh, "\n", sep=" ")
