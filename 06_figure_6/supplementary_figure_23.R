###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### LOAD LIBRARIES ---
library("stringr")
library("data.table")
library("reshape2")
library("pheatmap")
library("RColorBrewer")


### DEFINE PATH ---
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

dir.tfpeaks2gene <- file.path(dir.reproduce_data, "tfpeaks2gene")

### DEFINE FILES ---
file.rds_atacseq_masterdata <- file.path(dir.reproduce_data, "wcdt_atacseq_masterdata.rds")
file.rds_rnaseq <- file.path(dir.reproduce_data, "wcdt_rnaseq_atacseq_gene_expression.rds")
file.tf2genes <- file.path(dir.tfpeaks2gene, "wcdt_tf2genes_r_0p4_pval_0p05_distTSS_500kb.tsv.gz")
file.rds_genelist <- file.path(dir.reproduce_data, "genelist.rds")
file.utility_functions <- file.path(dir.reproduce_scripts, "00_prepare_data/00_01_utility_functions.R")


### LOAD FUNCTIONS ---
source(file.utility_functions)


### SAMPLEIDS ORDER AFTER CLUSTERING INDIVIDUAL SUBTYPES ---
sampleids <- c("DTB-127-BL","DTB-127-PRO","DTB-069-BL","DTB-172-BL","DTB-128-BL",
				"DTB-148-BL","DTB-090-PRO","DTB-167-BL","DTB-019-BL","DTB-019-PRO",
				"DTB-085-BL","DTB-111-PRO","DTB-149-BL","DTB-206-BL","DTB-129-BL",
				"DTB-176-BL","DTB-119-PRO","DTB-143-BL","DTB-101-BL","DTB-222-BL",
				"DTB-022-BL","DTB-156-BL","DTB-009-BL","DTB-035-BL","DTB-077-PRO",
				"DTB-146-BL","DTB-265-PRO","DTB-259-BL","PR-040-BL","DTB-092-BL",
				"PR-056-BL","DTB-055-PRO","PR-120-BL","DTB-135-BL","DTB-059-BL",
				"DTB-167-PRO","DTB-060-BL","DTB-138-BL","DTB-063-BL","DTB-080-BL",
				"DTB-126-BL","DTB-125-BL","DTB-011-BL","DTB-124-BL","DTB-140-BL",
				"DTB-141-BL","PR-095-BL","DTB-008-BL","DTB-074-BL","DTB-044-BL",
				"DTB-037-BL","DTB-091-BL","DTB-034-BL","DTB-004-BL","DTB-005-BL",
				"DTB-071-BL","DTB-064-BL","DTB-083-BL","DTB-053-BL","DTB-183-BL",
				"DTB-032-BL","PR-081-BL","DTB-036-BL","DTB-040-BL","DTB-061-BL",
				"DTB-080-PRO","DTB-218-BL","DTB-130-BL","DTB-097-BL","DTB-097-PRO")


### MOTIFS ---
motif_prefix_znf <- "MA0528.1_ZNF263"
motif_prefix_myc <- "MA0147.3_MYC"




############################################################################################################################################
### LOAD ATAC-SEQ MASTER DATA ---
list.rds_atacseq_masterdata <- readRDS(file=file.rds_atacseq_masterdata)

### GET METADATA ---
metadata <- list.rds_atacseq_masterdata$metadata_mcrpc
metadata <- metadata[match(sampleids, metadata$Sample_ID),]

### LOAD GENELIST ---
list.rds_genelist <- readRDS(file=file.rds_genelist)
genes_hallmark_myc <- unique( c( list.rds_genelist$msigdb_pathways$hallmark$HALLMARK_MYC_TARGETS_V1, list.rds_genelist$msigdb_pathways$hallmark$HALLMARK_MYC_TARGETS_V2) )



			
######################################################################
### TRIM USED DATA ---
des <- subset(metadata, select=c("Sample_ID","mCRPC_Subtype"))
des <- des %>% dplyr::mutate_all(as.character)
rownames(des) <- des$Sample_ID 
des$Sample_ID  <- NULL

### ITEMS LEVELS ---
items_subtypes <- c("AR+_NE-","ARL_NE-","AR+_NE+","AR-_NE+","AR-_NE-")

### COLOR ---
cp_subtypes <- c("#e31a1c","#fb9a99","#ff7f00","#33a02c","#1f78b4")

### ADD NAMES ---
names(cp_subtypes) <- items_subtypes

### DEFINE COLOR LIST ----
list.color <- list(mCRPC_Subtype=cp_subtypes)
######################################################################






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
######################################################################

## INVERSE NORMAL TRANSFORMATION ---
mat_norm <- as.matrix(t(apply(as.matrix(expr), 1, RankNorm)))





##########################################################################################################################
### LOAD TF2GENES DATA ---
dat_tf2genes <- data.table::fread(file=file.tf2genes, sep="\t", header=TRUE, nThread=50, data.table=FALSE, verbose=FALSE)

### GET MOTIF PEAKS ---
dat_znf <- dat_tf2genes[which(dat_tf2genes$motif_prefix == motif_prefix_znf),]

### SEPARATE CORR TARGET GENES ---
genes_znf_pos <- unique(dat_znf$Gene[which(dat_znf$Corr > 0)])

### ZNF263 TARGETS n MYC TARGETS GENES ---
genes_znf_hallmark_myc <- intersect(unique(dat_znf$Gene), genes_hallmark_myc)

##########################################################################################################################
### EXTRACT EXPRESSION DATA BY GENES ---
mat_norm_znf <- subset(mat_norm, rownames(mat_norm) %in% genes_znf_pos)
mat_norm_myc <- subset(mat_norm, rownames(mat_norm) %in% genes_znf_hallmark_myc)






##################################################################################################################
### GENERATE COLOR PALETTE ---
jColFun <- colorRampPalette(brewer.pal(n = 9, "RdBu"))

### HEATMAP: -----
file.plot1 <- file.path(dir.reproduce_fig, "supplementary_figure_23_a.pdf")
p1 <- pheatmap(mat_norm_znf, 
			color = rev(jColFun(128)), 
			kmeans_k = NA, breaks = NA, border_color = "#000000",
			cellwidth = NA, cellheight = NA, 
            scale = "none", 
			cluster_rows = TRUE, cluster_cols = TRUE, 
			clustering_distance_rows = "euclidean", 
			clustering_distance_cols = "euclidean", 
			clustering_method = "ward.D",
			cutree_rows = 2, cutree_cols = 3,
			treeheight_row = 25, treeheight_col = 25, 
			legend = TRUE, legend_breaks = NA, legend_labels = NA, 
			annotation_row = NA, annotation_col = des, 
			annotation_colors = list.color, annotation_legend = TRUE,
			annotation_names_row = FALSE, annotation_names_col = TRUE,
			drop_levels = TRUE, show_rownames = FALSE, show_colnames = FALSE, main = "",
			fontsize = 10, fontsize_row = 4, fontsize_col = 4,
			display_numbers = FALSE, number_format = "%.2f", number_color = "grey30", 
			fontsize_number = 0.8, 
            gaps_row = NULL, gaps_col = NULL, 
			labels_row = NULL, labels_col = NULL, 
            filename = file.plot1, width = 5, height = 3,
			silent = FALSE, na_col = "#DDDDDD")

### HEATMAP: -----
file.plot2 <- file.path(dir.reproduce_fig, "supplementary_figure_23_b.pdf")
p2 <- pheatmap(mat_norm_myc, 
			color = rev(jColFun(128)), 
			kmeans_k = NA, breaks = NA, border_color = "#000000",
			cellwidth = NA, cellheight = NA, 
            scale = "none", 
			cluster_rows = TRUE, cluster_cols = TRUE, 
			clustering_distance_rows = "euclidean", 
			clustering_distance_cols = "euclidean", 
			clustering_method = "ward.D",
			cutree_rows = 2, cutree_cols = 3,
			treeheight_row = 25, treeheight_col = 25, 
			legend = TRUE, legend_breaks = NA, legend_labels = NA, 
			annotation_row = NA, annotation_col = des, 
			annotation_colors = list.color, annotation_legend = TRUE,
			annotation_names_row = FALSE, annotation_names_col = TRUE,
			drop_levels = TRUE, show_rownames = FALSE, show_colnames = FALSE, main = "",
			fontsize = 10, fontsize_row = 4, fontsize_col = 4,
			display_numbers = FALSE, number_format = "%.2f", number_color = "grey30", 
			fontsize_number = 0.8, 
            gaps_row = NULL, gaps_col = NULL, 
			labels_row = NULL, labels_col = NULL, 
            filename = file.plot2, width = 5, height = 3,
			silent = FALSE, na_col = "#DDDDDD")

