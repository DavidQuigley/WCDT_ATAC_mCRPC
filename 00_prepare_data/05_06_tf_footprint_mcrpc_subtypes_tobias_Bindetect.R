###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### LOAD LIBRARIES ---
library("stringr")
library("data.table")
library("dplyr")

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
file.ref_genome <- file.path("/opt/reference/Homo_sapiens/NCBI/GRCh38Decoy/Sequence/WholeGenomeFasta/genome.fa") # HUMAN REF GENOME HG38
file.blacklist_region <- file.path(dir.reproduce, "reference/GRCh38_unified_blacklist_2020_05_05.bed") # ENCODE BLACKLIST HG38
file.motif_db <- file.path(dir.reproduce, "reference/JASPAR2018_CORE_vertebrates_non-redundant.meme")
file.header_peaks <- file.path(dir.bed, "peaks_header_annotation.txt")

### SUBTYPES ---
subtypes <- c("ARpNEn","ARpNEp","ARlNEn","ARnNEp","ARnNEn")

### GET PAIRWISE ALL GROUP COMBINATIONS ---
dm_grp <- do.call(rbind.data.frame, combn(x=subtypes, m=2, simplify=FALSE)) %>%
                dplyr::mutate_all(as.character)
colnames(dm_grp) <- c("Group1","Group2")
rownames(dm_grp) <- NULL

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
getBindetect <- function(dir.data, dir.output, file.motif_db, file.ref_genome, file.bed_peaks, file.header_peaks, group_id_1, group_id_2){
    # DEFINE FILES ---
    dir.output_group <- file.path(dir.output, paste("tobias_diff", paste(group_id_1, group_id_2, sep="_"), sep="_")  )
    dir.group_id_1 <- file.path(dir.data, group_id_1)
    dir.group_id_2 <- file.path(dir.data, group_id_2)
    file.bw_1 <- file.path(dir.group_id_1, paste(group_id_1, "footprints.bw", sep="_")) 
    file.bw_2 <- file.path(dir.group_id_2, paste(group_id_2, "footprints.bw", sep="_"))

    # DEFINE VARS ---
    ln_1 <-  paste("DIR_OUTPUT", dir.output_group, sep="=")
    ln_2 <-  paste("FILE_BIGWIG_1", file.bw_1, sep="=")
    ln_3 <-  paste("FILE_BIGWIG_2", file.bw_2, sep="=")
    ln_4 <-  paste("FILE_BED_PEAKS", file.bed_peaks, sep="=")
    ln_5 <-  paste("FILE_MOTIF_DB", file.motif_db, sep="=")
    ln_6 <-  paste("FILE_REF_GENOME", file.ref_genome, sep="=")
    ln_7 <-  paste("CONDITION_1", group_id_1, sep="=")
    ln_8 <-  paste("CONDITION_2", group_id_2, sep="=")
    ln_9 <-  paste("FILE_HEADER_PEAKS", file.header_peaks, sep="=")
    ln_10 <-  paste("DIR_OUTPUT", dir.output_group, sep="=")
    ln_10 <- "mkdir ${DIR_OUTPUT}"

    # TOBIAS FUNCTION ---
    line1 <- paste("TOBIAS BINDetect", "\\", sep=" ")
    line2 <- paste("  --signals ${FILE_BIGWIG_1} ${FILE_BIGWIG_2}", "\\", sep=" ")
    line3 <- paste("  --peaks ${FILE_BED_PEAKS}", "\\", sep=" ")
    line4 <- paste("  --peak-header ${FILE_HEADER_PEAKS}", "\\", sep=" ")
    line5 <- paste("  --motifs ${FILE_MOTIF_DB}", "\\", sep=" ")
    line6 <- paste("  --genome ${FILE_REF_GENOME}", "\\", sep=" ")
    line7 <- paste("  --cond-names ${CONDITION_1} ${CONDITION_2}", "\\", sep=" ")
    line8 <- paste("  --naming \"id_name\" ", "\\", sep=" ")
    line9 <- paste("  --motif-pvalue 1e-4", "\\", sep=" ")
    line10 <- paste("  --bound-pvalue 0.001", "\\", sep=" ")
    line11 <- paste("  --outdir ${DIR_OUTPUT}", "\\", sep=" ")
    line12 <- paste("  --cores 60", "\\", sep=" ")
    line14 <- paste("  --verbosity 5", sep=" ")

    # COMPILE COMMAND ---
    list.cmd <- list(ln_1, ln_2, ln_3, ln_4, ln_5, ln_6, ln_7, ln_8, ln_9, ln_10, 
                    line1, line2, line3, line4, line5, line6, line7, line8, line9, line10, line11, line12, line14)

    return(list.cmd)
}



### LOOP FOR EACH PAIR OF GROUPS ---
list.sh <- list()
for(k in 1:nrow(dm_grp)){
    group_id_1 <- dm_grp$Group1[k]
    group_id_2 <- dm_grp$Group2[k]

    id <- paste(group_id_1, group_id_2, sep="_")

    cat("START:", id, "\n", sep="\t")

    # GET BED FILE ---
    file.bed_peaks <- file.path(dir.bed, paste(id, ".bed", sep=""))

    # GET CMD ---
    list.cmd_sge <- get.SGEcmd(jobid=paste("tb", id, sep="_"), dir.logs)
    list.cmd_tobias <- getBindetect(dir.data=dir.footprints, dir.output=dir.footprints, file.motif_db, file.ref_genome, file.bed_peaks, file.header_peaks, group_id_1=group_id_1, group_id_2=group_id_2)

    # MERGE COMMANDS ---
    cmd <- c(unlist(list.cmd_sge), unlist(list.cmd_tobias))

    # WRITE QSUB SCRIPT --
    file.qsub <- file.path(dir.scripts, paste("qsub_tobias_bindetect_", id, ".sh", sep=""))
    write.table(cmd, file.qsub, row.names=FALSE, col.names=FALSE, quote=FALSE)

    list.sh[[k]] <- paste("qsub", file.qsub, sep=" ")

    cat("PROCESSED:", id, "\n", sep=" ")
}


### WRITE MAIN RUN SCRIPT ---
file.sh <- file.path(dir.scripts, "main_tobias_bindetect.sh")
write.table(unlist(list.sh), file.sh, row.names=FALSE, col.names=FALSE, quote=FALSE)

cat("FILE GENERATED:", file.sh, "\n", sep=" ")

