###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### DEFINE LIBRARIES ---
library("stringr")
library("data.table")
library("pheatmap")
library("RColorBrewer")

### DEFINE PATH ---
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

### DEFINE FILES ----
file.rds_atacseq_masterdata <- file.path(dir.reproduce_data, "wcdt_atacseq_masterdata.rds")
file.rds_atacseq_read_counts_mcrpc <- file.path(dir.reproduce_data, "atacseq_read_counts_mcrpc.rds")
file.kwtest <- file.path(dir.reproduce_data, "wcdt_atacseq_mcrpc_subtypes_KruskalWallisTest_summary.tsv.gz")
file.utility_functions <- file.path(dir.reproduce_scripts, "00_prepare_data/00_01_utility_functions.R")

### LOAD FUNCTIONS ---
source(file.utility_functions)


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


######################################################################
### LOAD ATAC-SEQ READ COUNTS DATA ---
list.rds_atacseq_read_counts_mcrpc <- readRDS(file=file.rds_atacseq_read_counts_mcrpc)
dat_atacseq <- list.rds_atacseq_read_counts_mcrpc$feature_counts_norm
rownames(dat_atacseq) <- dat_atacseq$feature_id
dat_atacseq$feature_id <- NULL
colnames(dat_atacseq) <- metadata$Sample_ID

### LOAD ATAC-SEQ PEAK ANNOTATION DATA ---
annot <- list.rds_atacseq_read_counts_mcrpc$feature_annotation


### LOAD KRUSKAL-WALLIS TEST DATA ---
dat_kwtest <- data.table::fread(file=file.kwtest, sep="\t", header=TRUE, nThread=50, data.table=FALSE, verbose=FALSE)
colnames(dat_kwtest)[1] <- "FeatureID"
dat_kwtest <- merge(dat_kwtest, annot, by="FeatureID")
dat_kwtest <- dat_kwtest[order(dat_kwtest$pvalue, decreasing=FALSE),]
dat_kwtest <- subset(dat_kwtest, dat_kwtest$pvalue <= 0.001)



### EXTRACT COUNTS DATA BY FEATUREIDS ---
mat <- subset(dat_atacseq, rownames(dat_atacseq) %in% dat_kwtest$FeatureID)

# NORMALIZE ---
mat <- as.matrix(mat)
mat[mat > 6] <- 6

### INVERSE NORMAL TRANSFORMATION ---
mat_norm <- as.matrix(t(apply(mat, 1, RankNorm)))


### SAMPLE ORDER FOR HEATMAP VISUALIZATION ----
sampleids <- c("DTB-009-BL","DTB-085-BL","DTB-148-BL","DTB-156-BL","DTB-077-PRO",
				"DTB-143-BL","DTB-149-BL","DTB-019-PRO","DTB-019-BL","DTB-022-BL",
				"DTB-111-PRO","DTB-172-BL","DTB-090-PRO","DTB-127-BL","DTB-127-PRO",
				"DTB-128-BL","DTB-101-BL","DTB-119-PRO","DTB-035-BL","DTB-146-BL",
				"DTB-222-BL","DTB-176-BL","DTB-129-BL","DTB-167-BL","DTB-069-BL",
				"DTB-206-BL","DTB-126-BL","DTB-064-BL","DTB-141-BL","DTB-091-BL",
				"DTB-037-BL","DTB-140-BL","DTB-125-BL","DTB-074-BL","DTB-092-BL",
				"DTB-044-BL","DTB-083-BL","DTB-011-BL","DTB-259-BL","DTB-008-BL",
				"DTB-124-BL","DTB-167-PRO","DTB-135-BL","DTB-005-BL","PR-056-BL",
				"DTB-059-BL","DTB-080-BL","DTB-063-BL","DTB-138-BL","DTB-034-BL",
				"PR-095-BL","DTB-060-BL","PR-120-BL","DTB-055-PRO","DTB-265-PRO",
				"PR-040-BL","DTB-071-BL","DTB-004-BL","DTB-053-BL","DTB-183-BL",
				"DTB-036-BL","PR-081-BL","DTB-040-BL","DTB-032-BL","DTB-218-BL",
				"DTB-097-BL","DTB-061-BL","DTB-080-PRO","DTB-130-BL","DTB-097-PRO")


### RE-ORDER SAMPLES ---
mat_norm <- subset(mat_norm, select=sampleids)


#### HEATMAP: -----
# GENERATE COLOR PALETTE ------
jColFun <- colorRampPalette(brewer.pal(n = 9, "RdBu"))

# PLOT ---
file.plot <- file.path(dir.reproduce_fig, "figure_3_a.pdf")
p <- pheatmap(mat_norm, 
			color = rev(jColFun(5)), 
			kmeans_k = NA, breaks = NA, border_color = NA,
			cellwidth = NA, cellheight = NA, 
            scale = "none", 
			cluster_rows = TRUE, cluster_cols = FALSE, 
			clustering_distance_rows = "euclidean", 
			clustering_distance_cols = "euclidean", 
			clustering_method = "ward.D2",
			cutree_rows = NA, cutree_cols = NA,
			treeheight_row = 25, treeheight_col = 25, 
			legend = TRUE, legend_breaks = NA, legend_labels = NA, 
			annotation_row = NA, annotation_col = des, 
			annotation_colors = list.color, annotation_legend = TRUE,
			annotation_names_row = FALSE, annotation_names_col = TRUE,
			drop_levels = TRUE, show_rownames = FALSE, show_colnames = TRUE, main = "",
			fontsize = 8, fontsize_row = 5, fontsize_col = 6,
			display_numbers = FALSE, number_format = "%.2f", number_color = "grey30", 
			fontsize_number = 0.8, 
            gaps_row = NULL, gaps_col = c(26,58,60,64), 
			labels_row = NULL, labels_col = NULL, 
            filename = file.plot, width = 10, height = 7,
			silent = FALSE, na_col = "#DDDDDD")



##### HEATMAP WITHOUT COLORBAR ---
#p <- pheatmap(mat_norm, 
#			color = rev(jColFun(5)),
#			kmeans_k = NA, breaks = NA, border_color = NA,
#			cellwidth = NA, cellheight = NA, 
#           scale = "none", 
#			cluster_rows = TRUE, cluster_cols = FALSE, 
#			clustering_distance_rows = "euclidean", 
#			clustering_distance_cols = "euclidean", 
#			clustering_method = "ward.D2",
#			cutree_rows = NA, cutree_cols = NA,
#			treeheight_row = 25, treeheight_col = 25, 
#			legend = FALSE, legend_breaks = NA, legend_labels = NA, 
#			annotation_row = NA, annotation_col = des, 
#			annotation_colors = list.color, annotation_legend = FALSE,
#			annotation_names_row = FALSE, annotation_names_col = TRUE,
#			drop_levels = TRUE, show_rownames = FALSE, show_colnames = FALSE, main = "",
#			fontsize = 8, fontsize_row = 5, fontsize_col = 5,
#			display_numbers = FALSE, number_format = "%.2f", number_color = "grey30", 
#			fontsize_number = 0.8, 
#           gaps_row = NULL, gaps_col = c(26,58,60,64), 
#			labels_row = NULL, labels_col = NULL, 
#            filename = file.plot, width = 3.5, height = 2,
#			silent = FALSE, na_col = "#DDDDDD")


