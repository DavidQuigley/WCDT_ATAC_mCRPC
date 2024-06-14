###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### DEFINE LIBRARIES ---
library("stringr")
library("data.table")

### DEFINE PATH ---
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/analysis/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

dir.footprints <- file.path(dir.reproduce_data, "tf_footprints/mcrpc_subtypes")

### DEFINE FILES ---
file.motifdb <- file.path(dir.reproduce, "reference/JASPAR2018_motif_gene_conversion.tsv")
file.rds_bindetect_results <- file.path(dir.footprints, "data_footprints/tobias_mcrpc_subtypes_bindetect_results.rds")

### SUBTYPES ---
subtypes <- c("ARpNEn","ARpNEp","ARlNEn","ARnNEp","ARnNEn")
subtypes_final <- c("AR+NE-","AR+NE+","ARlowNE-","AR-NE+","AR-NE-")




### GET PAIRWISE ALL GROUP COMBINATIONS ---
dm_grp <- do.call(rbind.data.frame, combn(x=subtypes, m=2, simplify=FALSE)) %>%
                dplyr::mutate_all(as.character)
colnames(dm_grp) <- c("Group1","Group2")
rownames(dm_grp) <- NULL
dm_grp$id <- apply(dm_grp, 1, function(x) paste(x[1], x[2], sep="_"))
dm_grp$f_name <- apply(dm_grp, 1, function(x) paste("tobias_diff", paste(x[1], x[2], sep="_"), sep="_"))



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

### FUNCTION: parseDataSummary() ---
parseDataSummary <- function(dat, motif_db, group_id_1, group_id_2){
    # ADD MOTIF PREFIX ---
    dat$motif_prefix <- NA
    for(i in 1:nrow(dat)){
        dat$motif_prefix[i] <- motif_db$motif_prefix[which(motif_db$motif_id == dat$motif_id[i])]
        dat$name[i] <- motif_db$motif_name[which(motif_db$motif_id == dat$motif_id[i])]
    }

    # GET COLUMNS OF INTEREST ---
    y_change <- colnames(dat)[stringr::str_detect(colnames(dat), "_change")]
    y_pvalue <- colnames(dat)[stringr::str_detect(colnames(dat), "_pvalue")]
    y_score <- colnames(dat)[stringr::str_detect(colnames(dat), "_mean_score")]
    y_bound <- colnames(dat)[stringr::str_detect(colnames(dat), "_bound")]

    # PREPARE DATA ---
    df <- subset(dat, select=c("motif_id","name","total_tfbs",y_score,y_bound,y_change,y_pvalue))
    colnames(df) <- c("motif_id","tf_name","total_tfbs",
                        "group1_mean_score","group2_mean_score",
                        "group1_bound","group2_bound",
                        "group1_group2_log2_foldchange",
                        "group1_group2_pvalue")


    # CORRECT PVALUE ---
    index0 <- which(df$group1_group2_pvalue == 0)
    if(length(index0) != 0){
        df$group1_group2_pvalue[index0] <- .Machine$double.eps #smallest value
    }

    # PVALUE CORRECTION ---
    df$fdr <- p.adjust(df$group1_group2_pvalue, method="fdr", n = length(df$group1_group2_pvalue))
    df$nlogfdr <- -log10(df$fdr)

    ### ADD GROUP INFO ---
    df <- cbind(group1=group_id_1, group2=group_id_2, df)
    rownames(df) <- NULL

    return(df)
}

### FUNCTION:getParams()  ---
getParams <- function(df, id){
    y_cutoff <- as.numeric(round(quantile(df$nlogfdr, 0.8), 0))

    score_threshold <- switch(id,  
                        "ARpNEn_ARpNEp" = 0.05,
                        "ARpNEn_ARlNEn" = 0.1,
                        "ARpNEn_ARnNEp" = 0.1,
                        "ARpNEn_ARnNEn" = 0.1,
                        "ARpNEp_ARlNEn" = 0.1,
                        "ARpNEp_ARnNEp" = 0.1,
                        "ARpNEp_ARnNEn" = 0.1,
                        "ARlNEn_ARnNEp" = 0.1,
                        "ARlNEn_ARnNEn" = 0.1,
                        "ARnNEp_ARnNEn" = 0.1)

    nlogp_threshold <- switch(id,  
                        "ARpNEn_ARpNEp" = y_cutoff,
                        "ARpNEn_ARlNEn" = y_cutoff,
                        "ARpNEn_ARnNEp" = y_cutoff,
                        "ARpNEn_ARnNEn" = y_cutoff,
                        "ARpNEp_ARlNEn" = y_cutoff,
                        "ARpNEp_ARnNEp" = y_cutoff,
                        "ARpNEp_ARnNEn" = y_cutoff,
                        "ARlNEn_ARnNEp" = y_cutoff,
                        "ARlNEn_ARnNEn" = y_cutoff,
                        "ARnNEp_ARnNEn" = y_cutoff)

    list.output <- list(score_threshold_pos=score_threshold,
                        score_threshold_neg=-score_threshold,
                        nlogp_threshold=nlogp_threshold)

    return(list.output)
}

### FUNCTION: filterData() ---
filterData <- function(df, id, group_id_1, group_id_2){
    # GET PARAMS ---
    list.params <- getParams(df, id)
    score_threshold_pos <- list.params$score_threshold_pos
    score_threshold_neg <- list.params$score_threshold_neg
    nlogp_threshold <- list.params$nlogp_threshold

    # TAG ANALYSIS ---
    #df$Analysis <- id

    # TAG SUBTYPE GROUP ---
    df$Group <- "GROUP_3"
    df$Group[which((df$group1_group2_log2_foldchange >= score_threshold_pos) & (df$nlogfdr >= nlogp_threshold))] <- group_id_1
    df$Group[which((df$group1_group2_log2_foldchange <= score_threshold_neg) & (df$nlogfdr >= nlogp_threshold))] <- group_id_2

    # FILDER DATA ---
    df_pass_1 <- subset(df, df$Group == group_id_1)
    df_pass_2 <- subset(df, df$Group == group_id_2)

    # COMBINE DATA ---
    dm <- rbind(df_pass_1, df_pass_2)
    rownames(dm) <- NULL

    dm$nlogfdr <- NULL
    dm$Group <- NULL

    return(dm)
}

### GET MOTIFS ---
motif_db <- ParseMotifDB(file.motifdb)

### LOAD BINDETECT RESULTS ---
list.bindetect_results <- readRDS(file=file.rds_bindetect_results)



### LOOP FOR EACH PAIR-WISE COMPARISON ---
list.dm <- list()
for(k in 1:nrow(dm_grp)){
    id <- dm_grp$id[k]
    f_name <- dm_grp$f_name[k]
    group_id_1 <- dm_grp$Group1[k]
    group_id_2 <- dm_grp$Group2[k]

    cat("STRAT:", id, "\n", sep=" ")  

    # PARSE DATA ---
    df <- parseDataSummary(dat=list.bindetect_results[[id]], motif_db, group_id_1, group_id_2)
    df$nlogfdr <- NULL
    list.dm[[id]] <- df

    # FILTER DATA ---
    #list.dm[[id]] <- filterData(df, id, group_id_1, group_id_2)
    
    cat("DONE:", id, "\n", sep=" ")  
}


### AGGREGATE DATA ---
dm <- do.call(rbind.data.frame, list.dm)
rownames(dm) <- NULL

### COLUMN NAMES ---
colnames(dm) <- c("group1","group2","motif_id","tf_name","total_tfbs",
					"group1_mean_score","group2_mean_score","group1_bound","group2_bound",
					"differential_binding_score_group1_group2_log2_foldchange",
					"differential_binding_pvalue",
					"differential_binding_qvalue_fdr")


### CHANGE mCRPC SUBTYPE ANNOTATION ---
for(i in 1:length(subtypes)){
    dm$group1[which(dm$group1 == subtypes[i])] <- subtypes_final[i]
    dm$group2[which(dm$group2 == subtypes[i])] <- subtypes_final[i]
}


### WRITE OUTPUT ---
file.tbl <- file.path(dir.reproduce_tbl, "supplementary_table_s4.tsv")
write.table(dm, file.tbl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

cat("FILE GENERATED:", file.tbl, "\n", sep=" ")

