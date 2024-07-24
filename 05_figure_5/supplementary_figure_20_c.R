###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### LOAD LIBRARIES ---
library("stringr")
library("data.table")
library("pheatmap")
library("RColorBrewer")

### DEFINE PATH ---
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

dir.tfpeaks2gene <- file.path(dir.reproduce_data, "tfpeaks2gene")

### DEFINE FILES ---
file.rds_atacseq_masterdata <- file.path(dir.reproduce_data, "wcdt_atacseq_masterdata.rds")
file.rds_atacseq_read_counts_mcrpc <- file.path(dir.reproduce_data, "atacseq_read_counts_mcrpc.rds")
file.rds_rnaseq <- file.path(dir.reproduce_data, "wcdt_rnaseq_atacseq_gene_expression.rds")
file.dat <- file.path(dir.tfpeaks2gene, "wcdt_tfpeaks2genes_correlated_r_0p4_pval_0p05_distTSS_500kb.tsv.gz")
file.utility_functions <- file.path(dir.reproduce_scripts, "00_prepare_data/00_01_utility_functions.R")

### LOAD FUNCTIONS ---
source(file.utility_functions)



############################################################################################################################################
### LOAD ATAC-SEQ MASTER DATA ---
list.rds_atacseq_masterdata <- readRDS(file=file.rds_atacseq_masterdata)

### GET METADATA ---
metadata <- list.rds_atacseq_masterdata$metadata_mcrpc



############################################################################################################################################
### LOAD ATAC-SEQ READ COUNTS DATA ---
list.rds_atacseq_read_counts_mcrpc <- readRDS(file=file.rds_atacseq_read_counts_mcrpc)

### GET NORMALIZED ATAC-SEQ READCOUNTS DATA ---
dat_atacseq <- list.rds_atacseq_read_counts_mcrpc$feature_counts_norm
rownames(dat_atacseq) <- dat_atacseq$feature_id
dat_atacseq$feature_id <- NULL

### INVERSE NORMAL TRANSFORMATION ---
dat_atacseq <- t(apply(dat_atacseq, 1, RankNorm))

### REPLACE BY SAMPLEIDS ---
colnames(dat_atacseq) <- metadata$Sample_ID



############################################################################################################################################
### LOAD RNA-SEQ GENE EXPRESSION OBJECT: TPM ---
cat(format(Sys.time(), "%b %d %X"), "LOADAING R DATA-OBJECT: RNA-seq EXPRESSION ... ", "\n", sep=" ")
    list.rnaseq <- readRDS(file=file.rds_rnaseq)

    # ANNOTATION ---
    gannot <- list.rnaseq$gene_annotation
    gannot <- subset(gannot, gannot$GeneType == "protein_coding")
    gannot <- gannot[-which(stringr::str_detect(gannot$EnsemblID, "PAR_Y")),]

    # EXPRESSION ---
    expr <- as.data.frame(list.rnaseq$gene_expression)
    expr <- subset(expr, expr$FEATURE_ID %in% gannot$EnsemblID)
    rownames(expr) <- gannot$Gene
    expr$FEATURE_ID <- NULL
    colnames(expr) <- stringr::str_replace(colnames(expr), "RP", "BL")

cat(format(Sys.time(), "%b %d %X"), "DONE! ", "\n", sep=" ")

### RE-ARRANGE EXPRESSION DATA ---
expr <- expr[order(rownames(expr), decreasing=FALSE),]
expr <- subset(expr, select=metadata$Sample_ID)

### LOG2 TRANSFORM [log2(TPM + 1)]---
expr <- log2(expr + 1)

### REMOVE GENES WITH STANDARD DEVIATION = 0 ---
std_dev <- apply(expr, 1, sd)
del.index_1 <- which(std_dev == 0)
expr <- expr[-del.index_1,]

### REMOVE GENES WITH VALUE TPM <= 0 IN > 25% OF SAMPLES ---
n <- round(ncol(expr) * (25/100), 0)
ln_0 <- apply(expr, 1, function(x) length(which(x <= log2(0+1) )))
del.index_2 <- which(ln_0 >= n)
expr <- expr[-del.index_2,]

### INVERSE NORMAL TRANSFORMATION ---
dat_rnaseq <- t(apply(expr, 1, RankNorm))



############################################################################################################################################
### LOAD TFPEAKS2GENE DATA ---
dat <- data.table::fread(file=file.dat, sep="\t", header=TRUE, nThread=60, data.table=FALSE, verbose=FALSE)
dat <- subset(dat, dat$Corr >= 0.5)
dat$FeatureID <- apply(dat, 1, function(x) paste(x[1], x[2], sep="_"))

### PREPARE ATAC-SEQ MATRIX ---
mat_atacseq <- matrix(NA, nrow=nrow(dat), ncol=ncol(dat_atacseq), dimnames=c(list(dat$FeatureID, colnames(dat_atacseq))))
mat_rnaseq <- matrix(NA, nrow=nrow(dat), ncol=ncol(dat_rnaseq), dimnames=c(list(dat$FeatureID, colnames(dat_rnaseq))))


### FILL MATRIX: ATACSEQ ---
n <- nrow(dat)
for(i in 1:n){
    featureid <- dat$Peak[i]
    index <- which(rownames(dat_atacseq) == featureid)
    mat_atacseq[i,] <- as.numeric(dat_atacseq[index,])
    cat("PROCESSED:", i, "OF", n, "\n", sep=" ")
}


### FILL MATRIX: RNASEQ ---
n <- nrow(dat)
for(i in 1:n){
    gene <- dat$Gene[i]
    index <- which(rownames(dat_rnaseq) == gene)
    mat_rnaseq[i,] <- as.numeric(dat_rnaseq[index,])
    cat("PROCESSED:", i, "OF", n, "\n", sep=" ")
}




############################################################################################################################################
######################################################################
### TRIM USED DATA ---
items <- c("Sample_ID","mCRPC_Subtype")

des <- subset(metadata, select=items)
des <- des %>% dplyr::mutate_all(as.character)

### MODIFY DESIGN TABLE ---
rownames(des) <- des$Sample_ID 

### ITEMS LEVELS ---
items_subtypes <- c("AR+_NE-","AR+_NE+","ARL_NE-","AR-_NE+","AR-_NE-")

### COLOR ---
cp_subtypes <- c("#e31a1c","#ff7f00","#fb9a99","#33a02c","#1f78b4")

### ADD NAMES ---
names(cp_subtypes) <- items_subtypes

### DEFINE COLOR LIST ----
list.color <- list(mCRPC_Subtype=cp_subtypes)
######################################################################




############################################################################################################################################
### FUNCTION: getClusterOrder() ---
getClusterOrder <- function(mat){
    obj_dist <- dist(x=mat, method="euclidean")
    obj_hclust <- hclust(d=obj_dist, method="ward.D2")
    y <- obj_hclust$labels[obj_hclust$order]
    return(y)
}

### GET CLUSTER ORDER ---
order_sampleid <- getClusterOrder(mat=t(mat_atacseq))
order_featureid <- getClusterOrder(mat=mat_atacseq)


### RE-ARRRANGE MATRIX ---
mat_atacseq <- mat_atacseq[match(order_featureid, rownames(mat_atacseq)),]
mat_atacseq <- mat_atacseq[,match(order_sampleid, colnames(mat_atacseq))]

mat_rnaseq <- mat_rnaseq[match(order_featureid, rownames(mat_rnaseq)),]
mat_rnaseq <- mat_rnaseq[,match(order_sampleid, colnames(mat_rnaseq))]


### RE-ARRANGE DESIGN TABLE ---
des <- des[match(order_sampleid, rownames(des)),]
des$Sample_ID  <- NULL



############################################################################################################################################
# GENERATE COLOR PALETTE ----
jColFun <- colorRampPalette(brewer.pal(n = 9, "RdBu"))

# PLOT: ATAC-SEQ ---
file.plot <- file.path(dir.reproduce_fig, "supplementary_figure_20_c_atacseq.pdf")
p1 <- pheatmap(mat_atacseq, 
			color = rev(jColFun(256)), 
			kmeans_k = NA, breaks = NA, border_color = NA,
			cellwidth = NA, cellheight = NA, 
            scale = "none", 
			cluster_rows = FALSE, cluster_cols = FALSE, 
			clustering_distance_rows = "euclidean", 
			clustering_distance_cols = "euclidean", 
			clustering_method = "ward.D2",
			cutree_rows = NA, cutree_cols = NA,
			treeheight_row = 25, treeheight_col = 25, 
			legend = TRUE, legend_breaks = NA, legend_labels = NA, 
			annotation_row = NA, annotation_col = des, 
			annotation_colors = list.color, annotation_legend = TRUE,
			annotation_names_row = FALSE, annotation_names_col = TRUE,
			drop_levels = TRUE, show_rownames = FALSE, show_colnames = FALSE, main = "",
			fontsize = 8, fontsize_row = 5, fontsize_col = 5,
			display_numbers = FALSE, number_format = "%.2f", number_color = "grey30", 
			fontsize_number = 0.8, 
            gaps_row = NULL, gaps_col = NULL, 
			labels_row = NULL, labels_col = NULL, 
            filename = file.plot, width = 4, height = 4,
			silent = FALSE, na_col = "#DDDDDD")


# PLOT: RNA-SEQ ---
file.plot <- file.path(dir.reproduce_fig, "supplementary_figure_20_c_rnaseq.pdf")
p2 <- pheatmap(mat_rnaseq, 
			color = rev(jColFun(256)), 
			kmeans_k = NA, breaks = NA, border_color = NA,
			cellwidth = NA, cellheight = NA, 
            scale = "none", 
			cluster_rows = FALSE, cluster_cols = FALSE, 
			clustering_distance_rows = "euclidean", 
			clustering_distance_cols = "euclidean", 
			clustering_method = "ward.D2",
			cutree_rows = NA, cutree_cols = NA,
			treeheight_row = 25, treeheight_col = 25, 
			legend = TRUE, legend_breaks = NA, legend_labels = NA, 
			annotation_row = NA, annotation_col = des, 
			annotation_colors = list.color, annotation_legend = TRUE,
			annotation_names_row = FALSE, annotation_names_col = TRUE,
			drop_levels = TRUE, show_rownames = FALSE, show_colnames = FALSE, main = "",
			fontsize = 8, fontsize_row = 5, fontsize_col = 5,
			display_numbers = FALSE, number_format = "%.2f", number_color = "grey30", 
			fontsize_number = 0.8, 
            gaps_row = NULL, gaps_col = NULL, 
			labels_row = NULL, labels_col = NULL, 
            filename = file.plot, width = 4, height = 4,
			silent = FALSE, na_col = "#DDDDDD")

