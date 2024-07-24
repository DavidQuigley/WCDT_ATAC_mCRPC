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
file.rds_atacseq_read_counts_tang2022_wcdt_combined <- file.path(dir.reproduce_data, "atacseq_read_counts_tang2022_wcdt_combined.rds")
file.rds_top_1percent_features_tang2022 <- file.path(dir.reproduce_data, "tang2022_atacseq_dqnorm_combat_top1percent.rds")
file.cluster_membership <- file.path(dir.reproduce_data, "wcdt_atacseq_unsupervised_cluster_membership.tsv") 
file.utility_functions <- file.path(dir.reproduce_scripts, "00_prepare_data/00_01_utility_functions.R")

### LOAD FUNCTIONS ---
source(file.utility_functions)



#####################################################################################################################

### LOAD ATAC-SEQ MASTER DATA ---
list.rds_atacseq_masterdata <- readRDS(file=file.rds_atacseq_masterdata)

### GET METADATA ---
metadata_tang2022_wcdt <- list.rds_atacseq_masterdata$metadata_tang2022_wcdt



### LOAD CLUSTER MEMBERSHIP DATA: FROM UNSUPERVISED CLUSTERING ---
dat.clust <- data.table::fread(file=file.cluster_membership, sep="\t", header=TRUE, nThread=1,  data.table=FALSE, verbose=FALSE)

### ADD CLUSTER INTO TO METADATA ---
metadata_tang2022_wcdt$Cluster <- NA
for(i in 1:nrow(dat.clust)){
    metadata_tang2022_wcdt$Cluster[ which( metadata_tang2022_wcdt$Sample_Name == dat.clust$SampleID[i] ) ] <- dat.clust$Cluster[i]
}

#####################################################################################################################


### LOAD FEATURES ---
mat_features <- readRDS(file=file.rds_top_1percent_features_tang2022)

### LOAD FEATURE COUNTS DATA ---
dat_dqnorm_combat <- readRDS(file=file.rds_atacseq_read_counts_tang2022_wcdt_combined)

### EXTRACT FEATURES ---
dat <- subset(dat_dqnorm_combat, rownames(dat_dqnorm_combat) %in% rownames(mat_features))
colnames(dat) <- metadata_tang2022_wcdt$Sample_Name_2


### GET SPEARMAN CORRELATION ---
mat_corr <- cor(dat, method="spearman")



######################################################################
### PREPARE DATA ---
des <- metadata_tang2022_wcdt %>% 
            dplyr::mutate_all(as.character)

### MODIFY DESIGN TABLE ---
des$Subtype_wcdt <- NA
des$Subtype_tang <- NA
des$Subtype_wcdt[which(des$Dataset == "WCDT")] <- des$Subtype[which(des$Dataset == "WCDT")]
des$Subtype_tang[which(des$Dataset == "TANG2022")] <- des$Subtype[which(des$Dataset == "TANG2022")]

des <- subset(des, select=c("Sample_Name_2","Subtype_tang","Subtype_wcdt","Cluster","Dataset"))

rownames(des) <- des$Sample_Name_2 
des$Sample_Name_2 <- NULL

### ITEMS LEVELS ---
#items_biopsy <- c("Bone","Lymph Node","Liver","Brain","CTC","Other")
items_subtypes_tang <- c("CRPC-AR","CRPC-NE","CRPC-SCL","CRPC-WNT")
items_subtypes_wcdt <- c("AR+_NE-","ARL_NE-","AR+_NE+","AR-_NE+","AR-_NE-")
items_cluster <- c("clust_1","clust_2","clust_3")
items_datasets <- c("WCDT","TANG2022")

### COLOR ---
#cp_biopsy <- c("#1f78b4","#e31a1c","#ff7f00","#525252","#ffff99","#cab2d6")
cp_subtypes_tang <- c("#e31a1c","#33a02c","#525252","#d9d9d9")
cp_subtypes_wcdt <- c("#e31a1c","#fb9a99","#ff7f00","#33a02c","#1f78b4")
cp_cluster <- c("#33a02c","#fb9a99","#e31a1c")
cp_datasets <- c("#e31a1c","#cab2d6")

### ADD NAMES ---
#names(cp_biopsy) <- items_biopsy
names(cp_subtypes_tang) <- items_subtypes_tang
names(cp_subtypes_wcdt) <- items_subtypes_wcdt
names(cp_cluster) <- items_cluster
names(cp_datasets) <- items_datasets

### DEFINE COLOR LIST ----
list.color <- list(Subtype_tang=cp_subtypes_tang,
					Subtype_wcdt=cp_subtypes_wcdt,
                    Cluster=cp_cluster,
					Dataset=cp_datasets)
#Biopsy_Site=cp_biopsy,
######################################################################



#################### HEATMAP #########################################

# BREAKLIST ---
breaksList <- seq(-0.2,1, by=0.1)

# GENERATE COLOR PALETTE ---
my_colors <- c("#dc2b1c","#ec6262","#f8dadd","#d2cde6","#5a58a3")
mypalette <- colorRampPalette(my_colors)(length(breaksList))



# PLOT ---
file.plot <- file.path(dir.reproduce_fig, "supplementary_figure_09.pdf")
p <- pheatmap(as.matrix(mat_corr), 
				color = rev(mypalette), 
				kmeans_k = NA, breaks = breaksList, border_color = NA,
				cellwidth = 1, cellheight = 1, 
            	scale = "none", 
				cluster_rows = TRUE, cluster_cols = TRUE, 
				clustering_distance_rows = "maximum", 
				clustering_distance_cols = "maximum", 
				clustering_method = "ward.D2",
				cutree_rows = 3, cutree_cols = 3,
				treeheight_row = 25, treeheight_col = 25, 
				legend = TRUE, legend_breaks = NA, legend_labels = NA, 
				annotation_row = NA, annotation_col = des, 
				annotation_colors = list.color, annotation_legend = TRUE,
				annotation_names_row = FALSE, annotation_names_col = TRUE,
				drop_levels = TRUE, show_rownames = FALSE, show_colnames = FALSE, main = "",
				fontsize = 6, fontsize_row = 3, fontsize_col = 3,
				display_numbers = FALSE, number_format = "%.2f", number_color = "grey30", 
				fontsize_number = 3, 
            	gaps_row = NULL, gaps_col = NULL, 
				labels_row = NULL, labels_col = NULL, 
            	filename = file.plot, width = 4, height = 4,
				silent = FALSE, na_col = "#DDDDDD")
