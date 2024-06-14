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

dir.analysis <- file.path(dir.wrk, "analysis/33_tf_nominate")
dir.data <- file.path(dir.analysis, "data")
dir.bigwig <- file.path(dir.wrk, "analysis/12_tf_footprint/subtypes_AR/tobias")
dir.bed <- file.path(dir.reproduce_data, "tf_footprints/mcrpc_subtypes/bed")
dir.output_bigwig <- file.path(dir.reproduce_data, "bigwig")

dir.temp <- tempdir() # REPLACE PATH

dir.scripts <- file.path(dir.reproduce, "scripts_temp")  # REPLACE PATH
dir.logs <- tempdir() # REPLACE PATH

### DEFINE FILES ---
file.chrom_sizes <- file.path(dir.wrk, "reference/grch38/grch38_sizes_major_chr.txt")
file.bw2wig <- file.path("/data1/home/rshrestha/software/utilities/bigWigToWig")
file.bedgraph2bigwig <- file.path("/data1/home/rshrestha/software/utilities/bedGraphToBigWig")

### SUBTYPES ---
subtype <- "ARnNEn"
tfs <- c("znf263","myc")


### FUNCTION: get.SGEcmd()---
get.SGEcmd <- function(jobid, dir.logs){
    jobid <- jobid
    line01 <- "#!/bin/bash"
    line02 <- "#$ -pe smp 12"
    line03 <- "#$ -V"
    line04 <- "#$ -R y"
    line05 <- "#$ -l mem_free=12G"
    line06 <- paste("#$ -N", jobid, sep=" ")
    line07 <- paste("#$ -o", file.path(dir.logs, paste(jobid, ".out", sep="")), sep=" ")
    line08 <- paste("#$ -e", file.path(dir.logs, paste(jobid, ".err", sep="")), sep=" ")

    # COMPILE CMD ---
    cmd <- list(line01, line02, line03, line04, line05, line06, line07, line08)

    return(cmd)    
}



### FUNCTION: getCMDs ---
getCMDs <- function(dir_bigwig, dir_footprint_bed, dir_output, dir_temp, file_chrom_sizes, file_bw2wig, file_bedgraph2bigwig, motif, subtype){
    # DEFINE CMDS ---
    #line_1 <- "#!/bin/bash"

    line_2 <- paste("MOTIF", motif, sep="=")
    line_3 <- paste("SUBTYPE", subtype, sep="=")

    line_4 <- paste("DIR_BIGWIG", dir_bigwig, sep="=")
    line_5 <- paste("DIR_FOOTPRINT_BED", dir_footprint_bed, sep="=")
    line_6 <- paste("DIR_OUTPUT", dir_output, sep="=")
    line_7 <- paste("DIR_TEMP", dir_temp, sep="=")

    line_8 <- paste("FILE_CHROM_SIZES", file_chrom_sizes, sep="=")
    line_9 <- paste("FILE_BW2WIG", file_bw2wig, sep="=")
    line_10 <- paste("FILE_BG2BW", file_bedgraph2bigwig, sep="=")

    line_11 <- "FILE_BIGWIG=${DIR_BIGWIG}/${SUBTYPE}/${SUBTYPE}_footprints.bw"
    line_12 <- "FILE_WIG=${DIR_TEMP}/${SUBTYPE}_tobias.wig"
    line_14 <- "FILE_WIG_BED=${DIR_TEMP}/${SUBTYPE}_tobias.bed"
    line_15 <- "FILE_BED_MOTIF=${DIR_FOOTPRINT_BED}/footprint_extended_${MOTIF}_${SUBTYPE}.bed"
    line_16 <- "FILE_WIG_BED_MOTIF=${DIR_TEMP}/footprint_subset_${MOTIF}_${SUBTYPE}.bed"
    line_17 <- "FILE_WIG_BED_MOTIF_SORTED=${DIR_TEMP}/footprint_subset_sorted_${MOTIF}_${SUBTYPE}.bed"
    line_18 <- "FILE_MOTIF_BEDGRAPH=${DIR_TEMP}/footprint_subset_${MOTIF}_${SUBTYPE}.bedgraph"
    line_19 <- "FILE_MOTIF_BIGWIG=${DIR_OUTPUT}/tobias_footprint_${MOTIF}_${SUBTYPE}.bw"

    line_20 <- "### 1. CONVERT BIGWIG TO WIG ---"
    line_21 <- "${FILE_BW2WIG} ${FILE_BIGWIG} ${FILE_WIG}"

    line_22 <- "### 2. CONVERT WIG TO BED ---"
    line_23 <- "wig2bed < ${FILE_WIG} > ${FILE_WIG_BED} --zero-indexed --max-mem=20G --sort-tmpdir=/scratch"

    line_24 <- "### 3. SUBSET WIG-BED BY FOOTPRINT BED ---"    
    line_25 <- "${BEDTOOLS}/intersectBed -a ${FILE_WIG_BED} -b ${FILE_BED_MOTIF} > ${FILE_WIG_BED_MOTIF} -u"

    line_26 <- "### 4. SORT BED ---"        
    line_27 <- paste("sort-bed --max-mem 20G --tmpdir",  file.path("/scratch"), "--unique ${FILE_WIG_BED_MOTIF} > ${FILE_WIG_BED_MOTIF_SORTED}", sep=" ")

    line_28 <- "### 5. CONVERT BED TO BEDGRAPH ---"        
    line_29 <- "bedops --element-of 1 ${FILE_WIG_BED} ${FILE_WIG_BED_MOTIF_SORTED} | awk \'BEGIN{ OFS=\"\\t\"; }{ print $1, $2, $3, $5 }\' > ${FILE_MOTIF_BEDGRAPH}"

    line_30 <- "### 6. CONVERT BEDGRAPH TO BIGWIG ---"        
    line_31 <- "${FILE_BG2BW} ${FILE_MOTIF_BEDGRAPH} ${FILE_CHROM_SIZES} ${FILE_MOTIF_BIGWIG}"

    # AGGREGATE CMDS ---
    list.cmd <- list(line_2, line_3, " ",  line_4, line_5, line_6, line_7, " ", line_8, line_9, line_10, " ", 
                     line_11, line_12, line_14, line_15, line_16, line_17, line_18, line_19, " ", 
                     line_20, line_21, " ",
                     line_22, line_23, " ",
                     line_24, line_25, " ",
                     line_26, line_27, " ",
                     line_28, line_29, " ",
                     line_30, line_31)
    
    return(list.cmd)                    
}


### PREPARE COMMANDS ---
list.sh <- list()
for(tf in tfs){
    key <- paste(tf, subtype, sep="_")

    # GET COMMANDS ---
    list.cmd_sge <- get.SGEcmd(jobid=paste("bw", key, sep="_"), dir.logs)
    list.cmd_bw <- getCMDs(dir_bigwig=dir.bigwig, 
                            dir_footprint_bed=dir.fp_bed, 
                            dir_output=dir.output_bigwig, 
                            dir_temp=dir.temp, 
                            file_chrom_sizes=file.chrom_sizes, 
                            file_bw2wig=file.bw2wig,
                            file_bedgraph2bigwig=file.bedgraph2bigwig, 
                            motif=tf, 
                            subtype=subtype)

    # MERGE COMMANDS ---
    cmd <- c(unlist(list.cmd_sge), unlist(list.cmd_bw))

    # WRITE QSUB SCRIPT --
    file.qsub <- file.path(dir.scripts, paste("qsub_", key, ".sh", sep=""))
    write.table(cmd, file.qsub, row.names=FALSE, col.names=FALSE, quote=FALSE)

    # STORE QSUB CMD ---
    list.sh[[key]] <- paste("qsub", file.qsub, sep=" ")

    cat("PROCESSED:", tf, subtype, "\n", sep=" ")
} 

### SUBMIT JOBS ---
for(i in 1:length(list.sh)){
    cmd <- list.sh[[i]]
    system(cmd)
}
