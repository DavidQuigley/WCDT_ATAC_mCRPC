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

#dir.footprints <- file.path(dir.wrk, "analysis/12_tf_footprint/subtypes_AR/tobias")
dir.footprints <- file.path(dir.reproduce_data, "tf_footprints/mcrpc_subtypes")
dir.bam <- file.path(dir.reproduce_data, "tf_footprints/mcrpc_subtypes/bam")
dir.bed <- file.path(dir.reproduce_data, "tf_footprints/mcrpc_subtypes/bed")

### DEFINE FILES ---
file.motifdb <- file.path(dir.reproduce, "reference/JASPAR2018_motif_gene_conversion.tsv")
file.rds_bindetect_results <- file.path(dir.footprints, "data_footprints/tobias_mcrpc_subtypes_bindetect_results.rds")

### SUBTYPES ---
subtypes <- c("ARpNEn","ARpNEp","ARlNEn","ARnNEp","ARnNEn")

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
parseDataSummary <- function(dat, motif_db){
    # ADD MOTIF PREFIX ---
    for(i in 1:nrow(dat)){
        dat$motif_prefix[i] <- motif_db$motif_prefix[which(motif_db$motif_id == dat$motif_id[i])]
        dat$name[i] <- motif_db$motif_name[which(motif_db$motif_id == dat$motif_id[i])]
    }

    # GET COLUMNS OF INTEREST ---
    y_change <- colnames(dat)[stringr::str_detect(colnames(dat), "_change")]
    y_pvalue <- colnames(dat)[stringr::str_detect(colnames(dat), "_pvalue")]

    # PREPARE DATA ---
    df <- subset(dat, select=c("motif_prefix","name","motif_id",y_change,y_pvalue))
    colnames(df) <- c("motif_prefix","motif_name","motif_id","score","pvalue")

    # CORRECT PVALUE ---
    index0 <- which(df$pvalue == 0)
    if(length(index0) != 0){
        df$pvalue[index0] <- .Machine$double.eps #smallest value
    }

    # PVALUE CORRECTION ---
    df$fdr <- p.adjust(df$pvalue, method="fdr", n = length(df$pvalue))
    df$nlogfdr <- -log10(df$fdr)

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
    df$Analysis <- id

    # TAG SUBTYPE GROUP ---
    df$Group <- "GROUP_3"
    df$Group[which((df$score >= score_threshold_pos) & (df$nlogfdr >= nlogp_threshold))] <- group_id_1
    df$Group[which((df$score <= score_threshold_neg) & (df$nlogfdr >= nlogp_threshold))] <- group_id_2

    # FILDER DATA ---
    df_pass_1 <- subset(df, df$Group == group_id_1)
    df_pass_2 <- subset(df, df$Group == group_id_2)

    # COMBINE DATA ---
    dm <- rbind(df_pass_1, df_pass_2)
    rownames(dm) <- NULL

    return(dm)
}


####################################################################################################################################
### GET PAIRWISE ALL GROUP COMBINATIONS ---
dm_grp <- do.call(rbind.data.frame, combn(x=subtypes, m=2, simplify=FALSE)) %>%
                dplyr::mutate_all(as.character)
colnames(dm_grp) <- c("Group1","Group2")
rownames(dm_grp) <- NULL
dm_grp$id <- apply(dm_grp, 1, function(x) paste(x[1], x[2], sep="_"))
dm_grp$f_name <- apply(dm_grp, 1, function(x) paste("tobias_diff", paste(x[1], x[2], sep="_"), sep="_"))



### LOAD BINDETECT RESULTS ----
list.bindetect_results <- readRDS(file=file.rds_bindetect_results)


### GET MOTIFS ---
motif_db <- ParseMotifDB(file.motifdb)


### LOOP FOR EACH PAIR-WISE COMPARISON ---
list.tobias_hits <- list()
for(k in 1:nrow(dm_grp)){
    id <- dm_grp$id[k]
    f_name <- dm_grp$f_name[k]
    subtype_1 <- dm_grp$Group1[k]
    subtype_2 <- dm_grp$Group2[k]

    cat("STRAT:", id, "\n", sep=" ")  

    # LOAD DATA ---
    df <- parseDataSummary(dat=list.bindetect_results[[id]], motif_db)

    # FILTER DATA ---
    list.tobias_hits[[id]] <- filterData(df, id, group_id_1=subtype_1, group_id_2=subtype_2)

    cat("DONE:", id, "\n", sep=" ")  
}

### AGGREGATE DATA ---
dat_tobias_hits <- do.call(rbind.data.frame, list.tobias_hits)
rownames(dat_tobias_hits) <- NULL

# -----
# dat_tobias_hits 
# -----



##################################################################################################################################
### COUNT MOTIF FREQUENCY BY SUBTYPE ---
dat_tobias_hits_freq <- dat_tobias_hits %>% 
                            dplyr::count(motif_prefix, Group, sort=FALSE, name="Freq")



### PREPARE MATRIX ---
mat_tobias_hits_freq <- matrix(0, nrow=nrow(motif_db), ncol=length(subtypes), dimnames=list(motif_db$motif_prefix, subtypes) )

### FILL MATRIX ---
for(i in 1:nrow(dat_tobias_hits_freq)){
    x <- dat_tobias_hits_freq$motif_prefix[i]
    y <- dat_tobias_hits_freq$Group[i]
    z <- dat_tobias_hits_freq$Freq[i]
    mat_tobias_hits_freq[x,y] <- z
}

### REMOVE MOTIFS WITHOUT HITS ---
mat_tobias_hits_freq <- mat_tobias_hits_freq[-which(rowSums(mat_tobias_hits_freq) == 0),]





### GET UNIQUE MOTIFSET ----
dat_tf_unq <- data.frame(motif_prefix=rownames(mat_tobias_hits_freq), gene=unlist(lapply(stringr::str_split(rownames(mat_tobias_hits_freq), "_"), function(x) x[2])))
genes.dup <- names(which(table(dat_tf_unq$gene) > 1))


### GET UNIQUESET OF MOTIFS/GENES ---
dat_tf_nodup <- subset(dat_tf_unq, !(dat_tf_unq$gene %in% genes.dup))
dat_tf_dup <- subset(dat_tf_unq, dat_tf_unq$gene %in% genes.dup)
dat_tf_dup <- dat_tf_dup[order(dat_tf_dup$gene),]

items <- c("ExtendedSite_AR","MA0476.1_FOS","MA1127.1_FOSB_JUN","MA0477.1_FOSL1","MA1130.1_FOSL2_JUN",
            "MA0655.1_JDP2","MA0488.1_JUN","MA0490.1_JUNB","MA0491.1_JUND","MA0056.1_MZF1","MA0161.2_NFIC",
            "MA0003.3_TFAP2A","MA0811.1_TFAP2B","MA0524.2_TFAP2C")
dat_tf_corrected <- subset(dat_tf_dup, dat_tf_dup$motif_prefix %in% items)

dat_tf <- rbind(dat_tf_nodup, dat_tf_corrected)
dat_tf <- dat_tf[order(dat_tf$gene),]
rownames(dat_tf) <- NULL


### EXTRACT UNIQUE SET OF MOTIF MATRIX ---
mat_tobias_hits_freq_unique <- subset(mat_tobias_hits_freq, rownames(mat_tobias_hits_freq) %in% dat_tf$motif_prefix)
dat_tobias_hits_freq_unique <- data.frame(cbind(motif_prefix=rownames(mat_tobias_hits_freq_unique), mat_tobias_hits_freq_unique))
dat_tobias_hits_freq_unique$gene <- unlist(lapply(stringr::str_split(dat_tobias_hits_freq_unique$motif_prefix, "_"), function(x) x[2]))
rownames(dat_tobias_hits_freq_unique) <- NULL




###### REFINE DATA AFTER CONFIRMATION OF TF EXPRESSION IN RNA-SEQ DATA PER SUBTYPE ------------------
file.mat_flagged <- file.path(dir.footprints, "data_footprints/tobias_mcrpc_subtypes_tf_footprints_hits_rnaseq_flagged.tsv")
mat_flagged <- data.table::fread(file=file.mat_flagged, sep="\t", header=TRUE, nThread=1, data.table=FALSE, verbose=FALSE)

### MOTIFS TO REMOVE AS THESE TF ARE NOT EXPRESSED IN RELAVANT SUBTYPE ---
ids_remove <- c("MA0914.1_ISL2","MA0710.1_NOTO","MA0890.1_GBX2","MA0667.1_MYF6","MA0845.1_FOXB1",
                "MA0717.1_RAX2","MA0618.1_LBX1","MA0820.1_FIGLA","MA0780.1_PAX3","MA0042.2_FOXI1",
                "MA0892.1_GSX1","MA0630.1_SHOX","MA0893.1_GSX2","MA0757.1_ONECUT3","MA0644.1_ESX1",
                "MA0718.1_RAX","MA0724.1_VENTX","MA0462.1_BATF_JUN","MA0726.1_VSX2","MA0680.1_PAX7",
                "MA0789.1_POU3F4","MA0723.1_VAX2")

### REMOVE TFs NOT EXPRESSED ---
mat_refine1 <- subset(dat_tobias_hits_freq_unique, !(dat_tobias_hits_freq_unique$motif_prefix %in% ids_remove))

### SEPARATE BY FLAGGED DATA ---
mat_refine2 <- rbind( subset(mat_refine1, !(mat_refine1$motif_prefix %in% mat_flagged$motif_prefix)), mat_flagged )





##################################################################################################################################
### COMPUTE METASCORE FOR EACH MOTIF-SUBTYPE -------
# meta_score = score * nlogfdr
dat_score <- dat_tobias_hits
dat_score$meta_score <- apply(dat_tobias_hits, 1, function(x) round(abs(as.numeric(x[4])) * as.numeric(x[7]), 4))

### ADD SCORE TO TF-HITS MATRIX
mat_new <- mat_refine2
rownames(mat_new) <- mat_new$motif_prefix
mat_new$motif_prefix <- NULL
mat_new$gene <- NULL


### CREATE EMPTY MATRIX ---
mat_agg <- matrix(0, nrow=nrow(mat_new), ncol=ncol(mat_new), dimnames=list(rownames(mat_new), colnames(mat_new)))

### GET AGGREGATE SCORE ---
for(i in 1:nrow(mat_new)){
    motif_prefix <- rownames(mat_new)[i]

    for(j in 1:ncol(mat_new)){
        subtype <- colnames(mat_new)[j]
        
        if(mat_new[i,j] != 0){
            d_temp <- dat_score[which((dat_score$motif_prefix == motif_prefix) & (dat_score$Group == subtype) ),]
            score_agg <- round( sum(d_temp$meta_score), 4)
            mat_agg[i,j] <- score_agg
        }

    }
}

### SCORE ---
mat_score <- data.frame( cbind(motif_prefix=rownames(mat_agg), mat_agg, gene=unlist(lapply(stringr::str_split(rownames(mat_agg), "_"), function(x) x[2]))) )
rownames(mat_score) <- NULL

mat_score$ARpNEn <- as.numeric(mat_score$ARpNEn)
mat_score$ARpNEp <- as.numeric(mat_score$ARpNEp)
mat_score$ARlNEn <- as.numeric(mat_score$ARlNEn)
mat_score$ARnNEp <- as.numeric(mat_score$ARnNEp)
mat_score$ARnNEn <- as.numeric(mat_score$ARnNEn)





###################################################################################################################################################
### ADD TO LIST ---
list.output <- list(tf_footprints_hits=dat_tobias_hits,
                    tf_footprints_hits_refined=mat_refine2,
                    tf_footprints_hits_score=mat_score)


### SAVE OBJECT TO RDATA FILE ---
file.rds_tf_hits <- file.path(dir.footprints, "data_footprints/tobias_mcrpc_subtypes_tf_hits.rds")
saveRDS(object=list.output, file=file.rds_tf_hits)
