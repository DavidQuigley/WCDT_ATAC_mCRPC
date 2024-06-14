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

### DEFINE FILE ---
file.motifdb <- file.path(dir.reproduce, "reference/JASPAR2018_motif_gene_conversion.tsv")
file.rds_atacseq_read_counts_mcrpc <- file.path(dir.reproduce_data, "atacseq_read_counts_mcrpc.rds")
file.tf_footprints_rds <- file.path(dir.footprints, "tobias_mcrpc_subtypes_tf_footprints_featureids_bound_unbound.rds")



###########################################################################################################################
### LOAD FOOTPRINT DATA ---
list.tf_footprints_rds <- readRDS(file=file.tf_footprints_rds)
list.footprints <- list.tf_footprints_rds$bound

### LOAD ATAC-SEQ READ COUNTS DATA ---
list.rds_atacseq_read_counts_mcrpc <- readRDS(file=file.rds_atacseq_read_counts_mcrpc)

### PREPARE ANNOTATION DATA ---
annot <- list.rds_atacseq_read_counts_mcrpc$feature_annotation
gr_annot <- GenomicRanges::makeGRangesFromDataFrame(df=annot, keep.extra.columns=TRUE)




###########################################################################################################################
### PREPARE MOTIF REGION GENOMIC RANGES ---
motifs_prefix <- names(list.footprints)
list.motif_pos <- list()
for(i in 1:length(motifs_prefix)){
    motif_prefix <- motifs_prefix[i]

    list.motif_pos_subtype <- list()
    # LOOP FOR EACH OF 5 AR-SUBTYPES ---
    for(j in 1:length( list.footprints[[motif_prefix]] )){
        subtype <- names(list.footprints[[motif_prefix]])[j]
        feature_ids <- list.footprints[[motif_prefix]][[subtype]]

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



###########################################################################################################################
### GET OVERLAP OF ATAC-SEQ PEAKS WITH MOTIF REGIONS ---
list.motif_peak <- list()
for(i in 1:length(list.motif_pos)){
    motif_prefix <- motifs_prefix[i]

    list.motif_peak_subtype <- list()
    # LOOP FOR EACH OF 5 AR-SUBTYPES ---
    for(j in 1:length( list.motif_pos[[motif_prefix]] )){
        subtype <- names( list.motif_pos[[motif_prefix]] )[j]
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

### WRITE OUTPUT ---
file.rds_tf2peaks <- file.path(dir.tfpeaks2gene, "wcdt_tf2peaks_mcrpc_subtypes.rds")
saveRDS(object=list.motif_peak, file=file.rds_tf2peaks)
