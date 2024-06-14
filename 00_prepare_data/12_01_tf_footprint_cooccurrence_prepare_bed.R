###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### CONFIGURANTION ---
options("scipen"=10)

### LOAD LIBRARIES ---
library("stringr")
library("data.table")
library("GenomicRanges")
library("plyranges")
library("dplyr")

### DEFINE PATH ---
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/analysis/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

dir.footprints <- file.path(dir.reproduce_data, "tf_footprints/mcrpc_subtypes")
dir.bed <- file.path(dir.reproduce_data, "tf_footprints/mcrpc_subtypes/bed")

### DEFINE FILES ---
file.rds_footprints <- file.path(dir.footprints, "data_footprints/tobias_mcrpc_subtypes_tf_footprints_featureids_bound_unbound.rds")

### SUBTYPES ---
subtype <- "ARnNEn"

### MOTIFS ---
motif_prefix_znf <- "MA0528.1_ZNF263"
motif_prefix_myc <- "MA0147.3_MYC"


##########################################################################################################################################################
### FUNCTION: getFootprintRanges() ---
getFootprintRanges <- function(feature_ids){
    # GET GENOMIC COORDINATES ---
    chr <- unlist(lapply(stringr::str_split(feature_ids, "_"), function(x) x[1]))
    start <- as.numeric(unlist(lapply(stringr::str_split(feature_ids, "_"), function(x) x[2])))
    end <- as.numeric(unlist(lapply(stringr::str_split(feature_ids, "_"), function(x) x[3])))
    
    # CREATE BED ---
    bed <- data.frame(chr=chr, start=start, end=end)

    # GET GENOMIC RANGES OBJECT ---
    gr <- GenomicRanges::makeGRangesFromDataFrame(bed, keep.extra.columns=TRUE)
    
    return(gr)
}

### FUNCTION: getExtendedRanges() ---
getExtendedRanges <- function(gr_motif){
    gr_extended <- dplyr::mutate(plyranges::anchor_center(gr_motif), width=6000)
    df_extended <- as.data.frame(gr_extended)
    bed <- subset(df_extended, select=c("seqnames","start","end"))

    return(bed)
}

##########################################################################################################################################################
### LOAD FOOTPRINTS ---
list.rds_footprints <- readRDS(file=file.rds_footprints)

### GET FOOTPRINTS BY MOTIF ---
footprint_featureids_znf <- list.rds_footprints$bound[[motif_prefix_znf]][[subtype]]
footprint_featureids_myc <- list.rds_footprints$bound[[motif_prefix_myc]][[subtype]]

### GET MOTIF FOOTPRINTS ---
gr_motif_znf <- getFootprintRanges(feature_ids=footprint_featureids_znf)
gr_motif_myc <- getFootprintRanges(feature_ids=footprint_featureids_myc)

### GET EXTENDED BED ---
bed_extended_znf <- getExtendedRanges(gr_motif=gr_motif_znf)
bed_extended_myc <- getExtendedRanges(gr_motif=gr_motif_myc)

### WRITE OUTPUT ---
file.bed_znf <- file.path(dir.bed, "footprint_extended_znf263_ARnNEn.bed")
file.bed_myc <- file.path(dir.bed, "footprint_extended_myc_ARnNEn.bed")

write.table(bed_extended_znf, file.bed_znf, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(bed_extended_myc, file.bed_myc, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

