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

### DEFINE FILES ----
file.rds_atacseq_masterdata <- file.path(dir.reproduce_data, "wcdt_atacseq_masterdata.rds")
file.utility_functions <- file.path(dir.reproduce_scripts, "00_prepare_data/00_01_utility_functions.R")

### LOAD FUNCTIONS ---
source(file.utility_functions)

### CHROMOSOMES ---
chromosome <- paste("chr", c(1:22,"X","Y"), sep="")

### GROUPS ---
subtypes <- c("ARpNEn","ARpNEp","ARlNEn","ARnNEp","ARnNEn")
group_ids <- c("AR+_NE-","AR+_NE+","ARL_NE-","AR-_NE+","AR-_NE-")



### FUNCTION: mergeBEDregions() ---
mergeBEDregions <- function(list.peaks_bed, sampleids){
    # GET NON-OVERLAPING REGIONS ---
    gr_nonoverlaping <- getNonOverlapPeaks(list.peaks=list.peaks_bed, sampleids=sampleids)
    gr_nonoverlaping <- sort(gr_nonoverlaping)


    # RETAIN PEAKS THAT ARE DETEDCTED IN AT LEAST 2 SAMPLES ---
    occurrences <- GenomicRanges::elementMetadata(gr_nonoverlaping) %>% 
                                        as.data.frame %>% 
                                        dplyr::select(-consensusIDs) %>% 
                                        rowSums

    gr_nonoverlaping <- gr_nonoverlaping[occurrences >= 2,]


    # GET DATA FRAME ---
    df.bed <- as.data.frame(gr_nonoverlaping)
    df.bed <- subset(df.bed, select=c("seqnames","start","end"))

    return(df.bed)
}


### FUNCTION: getBEDbyGroupPairs() ----
getBEDbyGroupPairs <- function(df_grp, subtype_1, subtype_2){
    items <- c("seqnames","start","end", subtype_1, subtype_2)
    d <- subset(df_grp, select=items)

    # REMOVE ROWS WITHOUT GROUP HITS ---
    del.index <- which( (d[,4] == 0) & (d[,5] == 0) )
    d <- d[-del.index,]

    # ADD KEY ---
    d$key <- apply(d, 1, function(x) paste(as.character(x[4]),as.character(x[5]), sep=":")) 

    # ADD STATUS ---
    d$status <- ifelse(d$key == "1:1", paste(subtype_1, subtype_2, sep=","), ifelse(d$key == "1:0", subtype_1, subtype_2))

    # EXTRACT DATA ---
    df_bed <- subset(d, select=c("seqnames","start","end","status"))
    rownames(df_bed) <- NULL

    return(df_bed)
}

### FUNCTION: getBEDbyGroup() ---
getBEDbyGroup <- function(df_nonoverlaping_regions, subtype){
    index_col <- which(colnames(df_nonoverlaping_regions) == subtype)
    index_row <- which(df_nonoverlaping_regions[,index_col] == 1)
    bed <- df_nonoverlaping_regions[index_row, 1:3]
    return(bed)
}






######################################################################
### LOAD ATAC-SEQ MASTER DATA ---
list.rds_atacseq_masterdata <- readRDS(file=file.rds_atacseq_masterdata)
metadata <- list.rds_atacseq_masterdata$metadata_mcrpc

### GET SUBTYPE SAMPLENAMES ---
list.samplenames <- list()
for(i in 1:length(group_ids) ){
    group_id <- group_ids[i]
    subtype <- subtypes[i]
    list.samplenames[[subtype]] <- metadata$Sample_Name[which(metadata$mCRPC_Subtype == group_id)]
}



######################################################################
### ATAC-SEQ PEAKS: UNCORRECTED --------
list.peaks <- list.rds_atacseq_masterdata$atacseq_peaks_uncorrected





######################################################################
### MERGE BED PER SUBTYPE ---
list.bed_subtype <- list()
for(subtype in subtypes){
    samplenames_subtype <- list.samplenames[[subtype]]
    list.peaks_subtype <- list.peaks[samplenames_subtype]

    # GET PEAK BED ---
    df <- mergeBEDregions(list.peaks_bed=list.peaks_subtype, sampleids=samplenames_subtype)

    # SORT ---
    df$seqnames <- as.character(df$seqnames)
    df <- df[order(df$seqnames, df$start),]
    rownames(df) <- NULL

    # ADD GROUP ID ---
    df$Subtype <- subtype

    # ADD TO LIST ---
    list.bed_subtype[[subtype]] <- df

    cat("DONE:", subtype, "\n", sep="\t")
}


### GET NON-OVERLAPING REGIONS ---
gr_nonoverlaping <- getNonOverlapPeaks(list.peaks=list.bed_subtype, sampleids=subtypes)
gr_nonoverlaping <- sort(gr_nonoverlaping)

### GET DATA FRAME ---
df_nonoverlaping_regions <- as.data.frame(gr_nonoverlaping)
df_nonoverlaping_regions$width <- NULL
df_nonoverlaping_regions$strand <- NULL


### WRITE OUTPUT: BED FORMAT ----
df.bed <- subset(df_nonoverlaping_regions, select=c("seqnames","start","end"))
file.output <- file.path(dir.reproduce_data, "tf_footprints/mcrpc_subtypes/bed/mcrpc_subtypes_atacseq_peaks_uncorrected_merged.bed")
write.table(df.bed, file.output, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)



######################################################################
### GET PAIRWISE GROUP COMBINATIONS ---
dm_grp <- do.call(rbind.data.frame, combn(x=subtypes, m=2, simplify=FALSE)) %>%
                dplyr::mutate_all(as.character)
colnames(dm_grp) <- c("Group1","Group2")
rownames(dm_grp) <- NULL


### LOOP FOR EACH PAIR OF GROUPS ---
for(k in 1:nrow(dm_grp)){
    subtype_1 <- dm_grp$Group1[k]
    subtype_2 <- dm_grp$Group2[k]

    id <- paste(subtype_1, subtype_2, sep="_")

    cat("START:", id, "\n", sep="\t")

    # GET BED BY PAIRS OF GROUPS ---
    df_bed <- getBEDbyGroupPairs(df_grp=df_nonoverlaping_regions, subtype_1, subtype_2)

    # WRITE OUTPUT ---
    file.output <- file.path(dir.reproduce_data, file.path("tf_footprints/mcrpc_subtypes/bed", paste(id, ".bed", sep="")  ))
    write.table(df_bed, file.output, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

    cat("DONE:", id, "\n", sep="\t")
}


### EXTRACT BED BY SUBTYPES ---
for(subtype in subtypes){
    bed <- getBEDbyGroup(df_nonoverlaping_regions, subtype)

    file.output <- file.path(dir.reproduce_data, "tf_footprints/mcrpc_subtypes/bed/mcrpc_subtypes_atacseq_peaks_uncorrected_subtype.bed")
    write.table(bed, file.output, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
}
