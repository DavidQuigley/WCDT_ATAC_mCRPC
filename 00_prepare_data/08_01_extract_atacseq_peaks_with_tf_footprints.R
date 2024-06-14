###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### DEFINE LIBRARIES ---
library("stringr")
library("data.table")
library("GenomicRanges")
library("IRanges")

### DEFINE PATH ---
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/analysis/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

dir.footprints <- file.path(dir.reproduce_data, "tf_footprints/mcrpc_subtypes")

### DEFINE FILES ---
file.motifdb <- file.path(dir.reproduce, "reference/JASPAR2018_motif_gene_conversion.tsv")
file.rds_atacseq_read_counts_mcrpc <- file.path(dir.reproduce_data, "atacseq_read_counts_mcrpc.rds")
file.rds_footprints <- file.path(dir.footprints, "data_footprints/tobias_mcrpc_subtypes_tf_footprints_featureids_bound_unbound.rds")

### SUBTYPES ---
subtypes <- c("ARpNEn","ARlNEn","ARpNEp","ARnNEp","ARnNEn")

### GET PAIRWISE ALL GROUP COMBINATIONS ---
dm_grp <- do.call(rbind.data.frame, combn(x=subtypes, m=2, simplify=FALSE)) %>%
                dplyr::mutate_all(as.character)
colnames(dm_grp) <- c("Group1","Group2")
rownames(dm_grp) <- NULL
dm_grp$id <- apply(dm_grp, 1, function(x) paste(x[1], x[2], sep="_"))



### FUNCTION: ParseMotifDB() ---
ParseMotifDB <- function(file.motifdb){
    # LOAD MOTIFDB DATA ---
    motifdb <- read.delim(file.motifdb, header=TRUE, stringsAsFactors=FALSE)
    motifdb <- subset(motifdb, motifdb$species == "Homo sapiens")
    
    motifdb$motif_prefix <- stringr::str_replace_all(motifdb$motif_prefix, "::", "")
    motifdb$motif_prefix <- stringr::str_replace_all(motifdb$motif_prefix, "[(]", "")
    motifdb$motif_prefix <- stringr::str_replace_all(motifdb$motif_prefix, "[)]", "")

    motif_prefix <- stringr::str_replace_all(motifdb$motif_name, "::", "_")
    motif_prefix <- stringr::str_replace_all(motif_prefix, "[(]", "_")
    motif_prefix <- stringr::str_replace_all(motif_prefix, "[)]", "")
    motif_prefix <- stringr::str_replace_all(motif_prefix, "[.]", "_")

    motifdb$motif_prefix_2 <- paste(motifdb$motif_id, motif_prefix, sep="_")
    motifdb <- motifdb[order(motifdb$motif_id, decreasing=FALSE),]

    rownames(motifdb) <- NULL

    motifdb <- subset(motifdb, select=c("motif_id","motif_name","motif_prefix_2"))
    colnames(motifdb) <- c("motif_id","motif_name","motif_prefix")
    motifdb <- motifdb[!duplicated(motifdb),]

    return(motifdb)
}



########################################################################################################################################
### GET MOTIFS ---
motif_db <- ParseMotifDB(file.motifdb)

### LOAD ATAC-SEQ PEAKS ANNOTATION ---
list.rds_atacseq_read_counts_mcrpc <- readRDS(file=file.rds_atacseq_read_counts_mcrpc)
annot <- list.rds_atacseq_read_counts_mcrpc$feature_annotation
annot <- subset(annot, select=c("seqnames","start","end","FeatureID"))
gr_annot <- GenomicRanges::makeGRangesFromDataFrame(df=annot, keep.extra.columns=TRUE)

### LOAD FOOTPRINTS ---
list.rds_footprints <- readRDS(file=file.rds_footprints)
motifs <- names(list.rds_footprints$bound)

### GET FOOTPRINT POSITION GRANGES ----
list.motif_pos <- list()
for(motif_prefix in motifs){
    list.motif_pos_subtype <- list()
    for(subtype in subtypes){
        feature_ids <- list.rds_footprints$bound[[motif_prefix]][[subtype]]

        list.feat <- stringr::str_split(feature_ids, "_")

        df_feat <- data.frame(chr=unlist(lapply(list.feat, function(x) x[1])),
                                start=as.numeric(unlist(lapply(list.feat, function(x) x[2]))),
                                end=as.numeric(unlist(lapply(list.feat, function(x) x[3]))),
                                motif_position=feature_ids,
                                motif_prefix=motif_prefix,
                                subtype=subtype)

        list.motif_pos_subtype[[subtype]] <- GenomicRanges::makeGRangesFromDataFrame(df=df_feat, keep.extra.columns=TRUE)
    }

    list.motif_pos[[motif_prefix]] <- list.motif_pos_subtype

    cat("PROCESSED:", motif_prefix, "\n", sep=" ")   
}




### GET OVERLAP OF ATAC-SEQ PEAKS WITH TF-FOOTPRINT REGIONS ---
list.motif_peak <- list()
for(motif_prefix in motifs){

    list.motif_peak_subtype <- list()
    # LOOP FOR EACH SUBTYPES ---
    for(subtype in subtypes){
        gr_motif <- list.motif_pos[[motif_prefix]][[subtype]]

        # GET OVERLAP OF ATAC-SEQ PEAK AND MOTIF REGION ---
        gr_overlap <- IRanges::mergeByOverlaps(query=gr_motif, subject=gr_annot)
        df_overlap <- as.data.frame(gr_overlap)

        items <- c("motif_prefix","subtype","motif_position","FeatureID")
        df_overlap <- subset(df_overlap, select=items)
        colnames(df_overlap) <- c("motif_prefix","subtype","motif_position","peak_position")

        list.motif_peak_subtype[[subtype]] <- df_overlap
    }

    list.motif_peak[[motif_prefix]] <- list.motif_peak_subtype
    cat("PROCESSED:", motif_prefix, "\n", sep=" ")  
}


### SAVE OBJECT TO RDATA FILE ---
file.rds_atacseq_peak_footprint_overlap <- file.path(dir.reproduce_data, "wcdt_atacseq_peak_footprint_overlap.rds")
saveRDS(object=list.motif_peak, file=file.rds_atacseq_peak_footprint_overlap)
