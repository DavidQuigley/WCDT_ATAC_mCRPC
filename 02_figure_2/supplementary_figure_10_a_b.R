###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### DEFINE LIBRARIES ---
library("stringr")
library("data.table")
library("dplyr")
library("pheatmap")
library("RColorBrewer")

### DEFINE PATH ---
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

### DEFINE FILES ----
file.rds_atacseq_masterdata <- file.path(dir.reproduce_data, "wcdt_atacseq_masterdata.rds")
file.rds_atacseq_read_counts_mcrpc <- file.path(dir.reproduce_data, "atacseq_read_counts_mcrpc.rds")
file.rds_rnaseq <- file.path(dir.reproduce_data, "wcdt_rnaseq_atacseq_gene_expression.rds")
file.rds_genelist <- file.path(dir.reproduce_data, "genelist.rds")
file.utility_functions <- file.path(dir.reproduce_scripts, "00_prepare_data/00_01_utility_functions.R")

# NOTE: THIS TUMOR CLUSTER DATA OBTAINED AFTER RUNNING UNSUPERVISED CLUSTERING OF ATAC-SEQ READCOUNT HEATMAP
# THIS DATA IS USED FOR LABELING THE CLUSTERS IN HEATMAP
file.cluster_membership <- file.path(dir.reproduce_data, "wcdt_atacseq_unsupervised_cluster_membership.tsv") 

### LOAD FUNCTIONS ---
source(file.utility_functions)



######################################################################
### GET GENELIST ---
list.rds_genelist <- readRDS(file=file.rds_genelist)
genes_labrecque2019 <- list.rds_genelist$genes_labrecque2019$Gene


######################################################################
### LOAD ATAC-SEQ MASTER DATA ---
list.rds_atacseq_masterdata <- readRDS(file=file.rds_atacseq_masterdata)
metadata <- list.rds_atacseq_masterdata$metadata_mcrpc

### EXTRACT DESIGN TABLE BY SUBTYPES ---
metadata1 <- metadata[which(metadata$mCRPC_Subtype == "AR+_NE-"),]
metadata2 <- metadata[which(metadata$mCRPC_Subtype == "ARL_NE-"),]
metadata3 <- metadata[which(metadata$mCRPC_Subtype == "AR+_NE+"),]
metadata4 <- metadata[which(metadata$mCRPC_Subtype == "AR-_NE+"),]
metadata5 <- metadata[which(metadata$mCRPC_Subtype == "AR-_NE-"),]

### COMBINE DESIGN TABLE ---
metadata <- rbind(metadata1, metadata2, metadata3, metadata4, metadata5)

######################################################################



######################################################################
### PREPARE DATA ---
des <- subset(metadata, select=c("Sample_ID","mCRPC_Subtype"))
des <- des %>% 
        dplyr::mutate_all(as.character)

### MODIFY DESIGN TABLE ---
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



################## RNA-SEQ DATA ################################################################################################

########## PREPARE RNA-SEQ EXPRESSION DATA -----------------------------------------------------------------------------------
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

########## -------------------------------------------------------------------------------------------------------------------


### NORMALIZE RNA-SEQ ---
expr[expr > 6] <- 6

### GET GENE INTERSECT ---
genes_rnaseq <- intersect(genes_labrecque2019, rownames(expr))

### EXTRACT EXPRESSION DATA BY GENES ---
expr <- subset(expr, rownames(expr) %in% genes_rnaseq)
expr <- expr[match(genes_rnaseq, rownames(expr)),]

### INVERSE NORMAL TRANSFORMATION ---
mat_expr_norm <- as.matrix(t(apply(expr, 1, RankNorm)))
mat_expr_norm <- as.matrix(mat_expr_norm)




#### HEATMAP: RNA-SEQ -----
# GENERATE COLOR PALETTE ------
jColFun <- colorRampPalette(brewer.pal(n = 9, "RdBu"))

# PLOT ---
file.plot1 <- file.path(dir.reproduce_fig, "supplementary_figure_10_a.pdf")
p1 <- pheatmap(mat_expr_norm, 
			color = rev(jColFun(256)), 
			kmeans_k = NA, breaks = NA, border_color = "#000000",
			cellwidth = 5, cellheight = 5, 
            scale = "none", 
			cluster_rows = FALSE, cluster_cols = FALSE, 
			clustering_distance_rows = "euclidean", 
			clustering_distance_cols = "euclidean", 
			clustering_method = "ward.D2",
			cutree_rows = 4, cutree_cols = 3,
			treeheight_row = 25, treeheight_col = 25, 
			legend = TRUE, legend_breaks = NA, legend_labels = NA, 
			annotation_row = NA, annotation_col = des, 
			annotation_colors = list.color, annotation_legend = TRUE,
			annotation_names_row = FALSE, annotation_names_col = TRUE,
			drop_levels = TRUE, show_rownames = TRUE, show_colnames = TRUE, main = "labrecque et al genes: RNA-seq",
			fontsize = 10, fontsize_row = 5, fontsize_col = 5,
			display_numbers = FALSE, number_format = "%.2f", number_color = "grey30", 
			fontsize_number = 0.8, 
            gaps_row = c(8,14,21), gaps_col = c(26,58,60,64), 
			labels_row = NULL, labels_col = NULL, 
            filename = file.plot1, width = 10, height = 4,
			silent = FALSE, na_col = "#DDDDDD")


################################################################################################################################



################## ATAC-SEQ DATA ################################################################################################

### LOAD ATAC-SEQ READ COUNTS DATA ---
list.rds_atacseq_read_counts_mcrpc <- readRDS(file=file.rds_atacseq_read_counts_mcrpc)
dat_atacseq <- list.rds_atacseq_read_counts_mcrpc$feature_counts_norm
dat_atacseq <- subset(dat_atacseq, select=c("feature_id", metadata$Sample_Name))
colnames(dat_atacseq) <- c("feature_id", metadata$Sample_ID)

### ATAC-SEQ PEAK ANNOTATION ---
annot_atacseq <- list.rds_atacseq_read_counts_mcrpc$feature_annotation
annot_atacseq <- subset(annot_atacseq, annot_atacseq$Gene %in% genes_labrecque2019)
annot_atacseq <- subset(annot_atacseq, annot_atacseq$Feature %in% c("Promoter","Intron","Distal Intergenic"))
annot_atacseq$Feature[which(annot_atacseq$Feature == "Distal Intergenic")] <- "Intergenic"
annot_atacseq$key <- apply(annot_atacseq, 1, function(x) paste(x[2], x[1], sep="_"))
annot_atacseq$GeneFeature <- apply(annot_atacseq, 1, function(x) paste(x[2], x[7], sep="_"))

### EXTRACT ATAC-SEQ READ COUNTS FEATURES BY GENES ---
dat_atacseq <- subset(dat_atacseq, dat_atacseq$feature_id %in% annot_atacseq$FeatureID)
dat_atacseq <- dat_atacseq[match(annot_atacseq$FeatureID, dat_atacseq$feature_id),]
dat_atacseq$feature_id <- annot_atacseq$GeneFeature





### FUNCTION: getunique_profile() ---
getunique_profile <- function(dat){
	unique.features <- sort(unique(as.character(dat[,1])), decreasing=FALSE)
	unique.features <- unique.features[unique.features != ""]
	dat <- subset(dat, dat[,1] != "")
	
	rown <- unique.features
	coln <- colnames(dat[2:ncol(dat)])
	mat <- matrix(nrow=length(unique.features), ncol=ncol(dat)-1, dimnames=list(rown, coln)) 

	ids <- as.character(dat[,1])
	for(i in 1:length(unique.features)){
		feature <- unique.features[i]
		index <- which(ids == feature)
		dat.sub <- as.matrix(dat[index,-1])
		
		if(nrow(dat.sub) == 1){
			mat[i,] <- as.numeric(dat.sub)
		} else {
			mat[i,] <- apply(dat.sub, 2, mean, na.rm=TRUE)
		}	
	}

	return(mat)
}


### GET UNIQUE ATAC-SEQ PROFILE ---
dat_atacseq_unique <- getunique_profile(dat_atacseq)



### ORDER ROWS BY FEATURES ----
features_ne1 <- c("CHGA_Intergenic","CHGA_Intron","CHGA_Promoter",
				"SYP_Promoter","ACTL6B_Intergenic","ACTL6B_Promoter",
				"SNAP25_Intron","SNAP25_Promoter",
				"INSM1_Intergenic","INSM1_Intron","INSM1_Promoter",
				"ASCL1_Promoter","CHRNB2_Promoter",
				"SRRM4_Intergenic","SRRM4_Intron","SRRM4_Promoter")

features_ne2 <- c("CELF3_Promoter","PCSK1_Intron","PCSK1_Promoter",
				"SOX2_Promoter","POU3F2_Intron",
				"POU3F2_Intergenic","POU3F2_Promoter",
				"LMO3_Intergenic","LMO3_Promoter",
				"NKX2-1_Intron","NKX2-1_Promoter")

features_ar <- c("AR_Intergenic","AR_Intron","AR_Promoter",
				"NKX3-1_Promoter","NKX3-1_Intergenic",
				"KLK3_Promoter","KLK3_Intron",
				"CHRNA2_Intron","CHRNA2_Promoter",
				"SLC45A3_Intergenic","SLC45A3_Intron","SLC45A3_Promoter",
				"TRGC1_Promoter","TRGC1_Intron",
				"NAP1L2_Intergenic","NAP1L2_Promoter",
				"S100A14_Promoter")

features_squam <- c("KRT6B_Intergenic","KRT6B_Promoter",
					"KRT6A_Intergenic","KRT6A_Promoter",
					"KRT5_Promoter","KRT5_Intergenic",
					"FGFBP1_Intergenic","FGFBP1_Promoter")
					
features_combine <- c(features_ne1, features_ne2, features_ar, features_squam)					
dat_atacseq_unique <- dat_atacseq_unique[match(features_combine, rownames(dat_atacseq_unique)),]


### NORMALIZE ---
dat_atacseq_unique[dat_atacseq_unique > 6] <- 6

### INVERSE NORMAL TRANSFORMATION ---
mat_norm_atacseq <- as.matrix(t(apply(dat_atacseq_unique, 1, RankNorm)))



#### HEATMAP: ATAC-SEQ -----
# GENERATE COLOR PALETTE ------
jColFun <- colorRampPalette(brewer.pal(n = 9, "RdBu"))

# PLOT ---
file.plot2 <- file.path(dir.reproduce_fig, "supplementary_figure_10_b.pdf")
p2 <- pheatmap(mat_norm_atacseq, 
			color = rev(jColFun(256)), 
			kmeans_k = NA, breaks = NA, border_color = "#000000",
			cellwidth = 5, cellheight = 5, 
            scale = "none", 
			cluster_rows = FALSE, cluster_cols = FALSE, 
			clustering_distance_rows = "euclidean", 
			clustering_distance_cols = "euclidean", 
			clustering_method = "ward.D2",
			cutree_rows = 4, cutree_cols = 3,
			treeheight_row = 25, treeheight_col = 25, 
			legend = TRUE, legend_breaks = NA, legend_labels = NA, 
			annotation_row = NA, annotation_col = des, 
			annotation_colors = list.color, annotation_legend = TRUE,
			annotation_names_row = FALSE, annotation_names_col = TRUE,
			drop_levels = TRUE, show_rownames = TRUE, show_colnames = TRUE, main = "labrecque et al genes: ATAC-seq",
			fontsize = 10, fontsize_row = 5, fontsize_col = 5,
			display_numbers = FALSE, number_format = "%.2f", number_color = "grey30", 
			fontsize_number = 0.8, 
            gaps_row = c(16,27,44), gaps_col = c(26,58,60,64), 
			labels_row = NULL, labels_col = NULL, 
            filename = file.plot2, width = 10, height = 6,
			silent = FALSE, na_col = "#DDDDDD")
