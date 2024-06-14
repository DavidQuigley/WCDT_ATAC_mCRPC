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
library("NbClust")
library("ggplot2")
library("gridExtra")


### DEFINE PATH ---
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/analysis/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

### DEFINE FILES ----
file.rds_atacseq_masterdata <- file.path(dir.reproduce_data, "wcdt_atacseq_masterdata.rds")
file.rds_atacseq_read_counts_mcrpc <- file.path(dir.reproduce_data, "atacseq_read_counts_mcrpc.rds")
file.utility_functions <- file.path(dir.reproduce_scripts, "00_prepare_data/00_01_utility_functions.R")

# NOTE: THIS TUMOR CLUSTER DATA OBTAINED AFTER RUNNING UNSUPERVISED CLUSTERING OF ATAC-SEQ READCOUNT HEATMAP
# THIS DATA IS USED FOR LABELING THE CLUSTERS IN HEATMAP
file.cluster_membership <- file.path(dir.reproduce_data, "wcdt_atacseq_unsupervised_cluster_membership.tsv") 

### LOAD FUNCTIONS ---
source(file.utility_functions)



######################################################################
### LOAD ATAC-SEQ MASTER DATA ---
list.rds_atacseq_masterdata <- readRDS(file=file.rds_atacseq_masterdata)

### LOAD ATAC-SEQ READ COUNTS DATA ---
list.rds_atacseq_read_counts_mcrpc <- readRDS(file=file.rds_atacseq_read_counts_mcrpc)


### GET METADATA ---
metadata <- list.rds_atacseq_masterdata$metadata_mcrpc

### GET WGS STATUS ---
status_wgs_mcrpc <- list.rds_atacseq_masterdata$status_wgs_mcrpc

### LOAD CLUSTER MEMBERSHIP DATA: FROM UNSUPERVISED CLUSTERING ---
dat.clust <- data.table::fread(file=file.cluster_membership, sep="\t", header=TRUE, nThread=1,  data.table=FALSE, verbose=FALSE)
dat.clust <- dat.clust[match(metadata$Sample_ID, dat.clust$SampleID),]
colnames(dat.clust) <- c("Sample_ID","Cluster")

### MERGE DATA ---
metadata <- merge(metadata, status_wgs_mcrpc, by.x="Sample_ID", by.y="SampleID")
metadata$Cluster <- dat.clust$Cluster


######################################################################
#################### PREPARE HEADER DATA #############################
### TRIM USED DATA ---
items <- c("Sample_ID","Cluster","mCRPC_Subtype","Biopsy_Site",
			"CDK12_2loss","CHD1_2loss","SPOP_inact_miss",
			"RB1_2loss","PTEN_2loss",
			"TMPRSS2_ERG_fusion","MYC_amp",
			"AR_amp","AR_enhancer_amp")

des <- subset(metadata, select=items)
des[is.na(des)] <- FALSE

### ADJUST FOR SAMPLE WITH MISSING WGS DATA ---------
id <- "DTB-004-BL"
des$AR_amp[which(des$Sample_ID == id)] <- NA
des$AR_enhancer_amp[which(des$Sample_ID == id)] <- NA
des$MYC_amp[which(des$Sample_ID == id)] <- NA
des$TMPRSS2_ERG_fusion[which(des$Sample_ID == id)] <- NA 
des$PTEN_2loss[which(des$Sample_ID == id)] <- NA
des$RB1_2loss[which(des$Sample_ID == id)] <- NA
des$SPOP_inact_miss[(des$Sample_ID == id)] <- NA
des$CHD1_2loss[which(des$Sample_ID == id)] <- NA
des$CDK12_2loss[which(des$Sample_ID == id)] <- NA


### PREPARE HEADER DATA ---
des <- des %>% dplyr::mutate_all(as.character)
rownames(des) <- des$Sample_ID 
des$Sample_ID  <- NULL




######################################################################
#################### PREPARE HEADER ANNOTATION #######################
### ITEMS LEVELS ---
items_cluster <- c("clust_1","clust_2","clust_3")
items_subtypes_mcrpc <- c("AR+_NE-","ARL_NE-","AR+_NE+","AR-_NE+","AR-_NE-")
items_biopsy_site <- c("Bone","Liver","Lymph Node","Other")
items_mut_cdk12 <- c("TRUE","FALSE")
items_loss_chd1 <- c("TRUE","FALSE")
items_mut_spop <- c("TRUE","FALSE")
items_loss_rb1 <- c("TRUE","FALSE")
items_loss_pten <- c("TRUE","FALSE")
items_fusion_t2e <- c("TRUE","FALSE")
items_amp_myc <- c("TRUE","FALSE")
items_amp_ar <- c("TRUE","FALSE")
items_amp_ar_enh <- c("TRUE","FALSE")


### COLOR ---
cp_cluster <- c("#33a02c","#fb9a99","#e31a1c")
cp_subtypes_mcrpc <- c("#e31a1c","#fb9a99","#ff7f00","#33a02c","#1f78b4")
cp_biopsy_site <- c("#1f78b4","#ff7f00","#e31a1c","#cab2d6")
cp_mut_cdk12 <- c("#33a02c","#bdbdbd")
cp_loss_chd1 <- c("#1f78b4","#bdbdbd")
cp_mut_spop <- c("#33a02c","#bdbdbd")
cp_loss_rb1 <- c("#1f78b4","#bdbdbd")
cp_loss_pten <- c("#1f78b4","#bdbdbd")
cp_fusion_t2e <- c("#000000","#bdbdbd")
cp_amp_myc <- c("#e31a1c","#bdbdbd")
cp_amp_ar <- c("#e31a1c","#bdbdbd")
cp_amp_ar_enh <- c("#e31a1c","#bdbdbd")

### ADD NAMES ---
names(cp_cluster) <- items_cluster
names(cp_subtypes_mcrpc) <- items_subtypes_mcrpc
names(cp_biopsy_site) <- items_biopsy_site
names(cp_mut_cdk12) <- items_mut_cdk12
names(cp_loss_chd1) <- items_loss_chd1
names(cp_mut_spop) <- items_mut_spop
names(cp_loss_rb1) <- items_loss_rb1
names(cp_loss_pten) <- items_loss_pten
names(cp_fusion_t2e) <- items_fusion_t2e
names(cp_amp_myc) <- items_amp_myc
names(cp_amp_ar) <- items_amp_ar
names(cp_amp_ar_enh) <- items_amp_ar_enh

### DEFINE COLOR LIST ----
list.color <- list(Cluster=cp_cluster,
					mCRPC_Subtype=cp_subtypes_mcrpc,
					Biopsy_Site=cp_biopsy_site,
					CDK12_2loss=cp_mut_cdk12,
					CHD1_2loss=cp_loss_chd1,
					SPOP_inact_miss=cp_mut_spop,
					RB1_2loss=cp_loss_rb1,
					PTEN_2loss=cp_loss_pten,
					TMPRSS2_ERG_fusion=cp_fusion_t2e,
					MYC_amp=cp_amp_myc,
					AR_amp=cp_amp_ar,
					AR_enhancer_amp=cp_amp_ar_enh)

######################################################################







######################################################################
#################### PREPARE ATAC-SEQ MATRIX #########################

### GET NORMALIZED ATAC-SEQ READCOUNTS DATA ---
dat <- list.rds_atacseq_read_counts_mcrpc$feature_counts_norm
rownames(dat) <- dat$feature_id
dat$feature_id <- NULL

#> all(colnames(dat) == metadata$Sample_Name)
#[1] TRUE

colnames(dat) <- metadata$Sample_ID


### EXTRACT FEATURES --- # LOADED FROM FILE.utility_functions #
# 25% most variable ATAC-seq peaks
list.features_filter <- featurefilter(mydata=dat, percentile=25)
mat <- list.features_filter$filtered_data

### GET SPEARMAN CORRELATION ---
mat_corr <- cor(mat, method="spearman")

######################################################################





######################################################################
#################### HEATMAP #########################################

# BREAKLIST ---
breaksList <- seq(0.3,1, by=0.1)

# GENERATE COLOR PALETTE ---
my_colors <- c("#dc2b1c","#ec6262","#f8dadd","#d2cde6","#5a58a3")
mypalette <- colorRampPalette(my_colors)(length(breaksList))



# PLOT ---
file.plot <- file.path(dir.reproduce_fig, "figure_2_a.pdf")
p <- pheatmap(as.matrix(mat_corr), 
				color = rev(mypalette), 
				kmeans_k = NA, breaks = breaksList, border_color = NA,
				cellwidth = 2, cellheight = 2, 
            	scale = "none", 
				cluster_rows = TRUE, cluster_cols = TRUE, 
				clustering_distance_rows = "euclidean", 
				clustering_distance_cols = "euclidean", 
				clustering_method = "ward.D2",
				cutree_rows = 3, cutree_cols = 3,
				treeheight_row = 5, treeheight_col = 5, 
				legend = TRUE, legend_breaks = NA, legend_labels = NA, 
				annotation_row = NA, annotation_col = des, 
				annotation_colors = list.color, annotation_legend = FALSE,
				annotation_names_row = FALSE, annotation_names_col = TRUE,
				drop_levels = TRUE, show_rownames = TRUE, show_colnames = TRUE, main = "",
				fontsize = 3, fontsize_row = 2, fontsize_col = 2,
				display_numbers = FALSE, number_format = "%.2f", number_color = "grey30", 
				fontsize_number = 3, 
            	gaps_row = NULL, gaps_col = NULL, 
				labels_row = NULL, labels_col = NULL, 
            	filename = file.plot, width = 4, height = 4,
				silent = FALSE, na_col = "#DDDDDD")

######################################################################




######################################################################
#################### PARSE SAMPLE CLUSTERS ###########################
### WRITE CLUSTER OUTPUT ---
sampleids_order <- p$tree_col$labels[p$tree_col$order]
cluster <- c(rep("clust_1",8), rep("clust_2",24), rep("clust_3",38))
df_clust <- data.frame( cbind(SampleID=sampleids_order, Cluster=cluster) )

file.cluster_membership <- file.path(dir.reproduce_data, "wcdt_atacseq_unsupervised_cluster_membership.tsv") 
#write.table(df_clust, file.cluster_membership, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)






######### SUPPLEMENTARY FIGURE S6 #######################################################################
######### ESTIMATE OPTIMUM NUMBER OF CLUSTERS USING ELBOW METHOD ######################################## 
### PCA ---
obj_pca <- prcomp(mat_corr)

### ESTIMATE OPTIMUM NUMBER CLUSTERS  -----
set.seed(123)
res_nbclust <- NbClust::NbClust(data = obj_pca$rotation[,1:3], diss = NULL, distance = "euclidean", 
                    		min.nc = 2, max.nc = 10, 
                    		method = "ward.D2", index = "all", 
                    		alphaBeale = 0.1)

### PLOT: ELBOW METHOD ---
p_clust <- factoextra::fviz_nbclust(x=obj_pca$rotation[,1:3], FUNcluster=kmeans, method="wss")  +
                        geom_vline(xintercept = 3, linetype = 2, color="steelblue") +
                        #labs(subtitle = "Elbow method", size=6) +
            			theme(
                			axis.text.x = element_text(size = 6, color="#000000"),
                			axis.text.y = element_text(size = 6, color="#000000"),
                			axis.title = element_text(size = 6, color="#000000"),
                			plot.title = element_text(size = 6, color="#000000", hjust=0.5),	
							plot.subtitle = element_text(size = 6, color="#000000", hjust=0.5),
                			panel.grid.major = element_blank(),
                			panel.grid.minor = element_blank(),
                			axis.ticks = element_line(size=0.2, color="#000000"),	
                			strip.text = element_text(size=6, color="#000000"),
                			strip.background = element_rect(fill="#FFFFFF", color="#FFFFFF"),
                			panel.background = element_rect(fill="#FFFFFF", color="#000000"),
                			legend.text = element_text(size = 6, color="#000000"),
                			legend.title = element_blank(),
                			legend.key.size = unit(0.3, "cm"),
                			legend.position = "none", legend.box = "horizontal") +
            			guides(fill = guide_legend(nrow = 1)) + 
            			ylab("Total Within Sum of Square") +
            			xlab("Number of clusters k") + 
						ggtitle("") 
            			#ggtitle("Optimal number of clusters", subtitle = "Elbow method") 


### GENERATE PLOT ---
file.plot <- file.path(dir.reproduce_fig, "supplementary_figure_08.pdf")
pdf(file.plot, width=2, height=2)
    grid.arrange(p_clust, nrow=1, ncol=1)
dev.off()





########### FISHER'S EXACT TEST: ATAC-seq Clusters vs mCRPC clusters ###########
df <- subset(metadata, select=c("Cluster","mCRPC_Subtype"))

#> table(df$Cluster, df$mCRPC_Subtype)
#
#          AR-_NE- AR-_NE+ AR+_NE- AR+_NE+ ARL_NE-
#  clust_1       3       3       0       0       2
#  clust_2       3       0       6       1      14
#  clust_3       0       1      20       1      16


### FUNCTION: getContingencyMatrix() --
getContingencyMatrix <- function(df, cluster, subtype){
	dm <- df 
	dm$Cluster <- ifelse(df$Cluster == cluster, 1, 0)
	dm$mCRPC_Subtype <- ifelse(df$mCRPC_Subtype == subtype, 1, 0)

	d <- dm %>% 
			dplyr::group_by(Cluster, mCRPC_Subtype) %>% 
			dplyr::tally()

	mat <- matrix(rev(d$n), nrow=2, ncol=2, byrow=TRUE, dimnames=list(c(1,0),c(1,0)))

	return(mat)
}



### CLUSTER-3 vs AR+NE-
mat <- getContingencyMatrix(df, cluster="clust_3", subtype="AR+_NE-")
fisher.test(x=mat, alternative="two.sided")

#> fisher.test(x=mat, alternative="two.sided")
#
#        Fisher's Exact Test for Count Data
#
#data:  mat
#p-value = 0.005839
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.456312 17.300565
#sample estimates:
#odds ratio
#   4.70256



### CLUSTER-2 vs ARlNE-  ------
mat <- getContingencyMatrix(df, cluster="clust_2", subtype="ARL_NE-")
fisher.test(x=mat, alternative="two.sided")

#> fisher.test(getContingencyMatrix(df, cluster="clust_2", subtype="ARL_NE-"))
#
#        Fisher's Exact Test for Count Data
#
#data:  getContingencyMatrix(df, cluster = "clust_2", subtype = "ARL_NE-")
#p-value = 0.1393
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
# 0.7129219 6.7398807
#sample estimates:
#odds ratio
#  2.153024

### CLUSTER-1 vs AR-NE+ ------
mat <- getContingencyMatrix(df, cluster="clust_1", subtype="AR-_NE+")
fisher.test(x=mat, alternative="two.sided")

#> fisher.test(getContingencyMatrix(df, cluster="clust_1", subtype="AR-_NE+"))
#
#        Fisher's Exact Test for Count Data
#
#data:  getContingencyMatrix(df, cluster = "clust_1", subtype = "AR-_NE+")
#p-value = 0.003863
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#    2.182988 1919.419979
#sample estimates:
#odds ratio
#   32.3312

#############################################################################################################

###### FISHERS' EXACT TEST --------------------------------------
items_clusters <- sort(unique(metadata$Cluster))
items_biopsy_sites <- sort(unique(metadata$Biopsy_Site))

ctr <- 1
list.dm <- list()
for(clust in items_clusters){
	for(biopsy_site in items_biopsy_sites){
		df <- data.frame(SampleID=metadata$Sample_ID, Cluster=0, BiopsySite=0)
		df$Cluster[ which(metadata$Cluster == clust) ] <- 1
		df$BiopsySite[ which(metadata$Biopsy_Site == biopsy_site) ] <- 1

		# GET TABULAR MATRIX ---
		tbl <- as.matrix(table(df$Cluster, df$BiopsySite))
		tbl1 <- c(tbl[2,2], tbl[2,1], tbl[1,2], tbl[1,1])

		# CONFUSION MATRIX ---
		mat_conf <- matrix( tbl1, byrow=TRUE, nrow=2, ncol=2, dimnames=list(c("clust_yes","clust_no"),c("bsite_yes","bsite_no")))

		# FISHERS' EXACT TEST ---
		res_fishers_test <- fisher.test(x=mat_conf)

		list.dm[[ctr]] <- data.frame(Cluster=clust, 
										BiopsySite=biopsy_site, 
										Clust_yes_BSite_yes=tbl[2,2], 
										Clust_yes_BSite_no=tbl[2,1], 
										Clust_no_BSite_yes=tbl[1,2], 
										Clust_no_BSite_no=tbl[1,1],
										pvalue=res_fishers_test$p.value)
		ctr <- ctr + 1
	}
}

### AGGREGATE DATA ---
dm <- do.call(rbind.data.frame, list.dm)
rownames(dm) <- NULL

#> subset(dm, dm$pvalue <= 0.05)
#   Cluster BiopsySite Clust_yes_BSite_yes Clust_yes_BSite_no Clust_no_BSite_yes Clust_no_BSite_no       pvalue
#1  clust_1       Bone                   0                  8                 28                34 0.0180997477
#5  clust_2       Bone                  17                  7                 11                35 0.0002455397 #ASSOCIATION
#7  clust_2 Lymph Node                   2                 22                 22                24 0.0011501815
#11 clust_3 Lymph Node                  18                 20                  6                26 0.0219279605





### FUNCTION: getFishersExactTest() ---
getFishersExactTest <- function(dat){
	items_clusters <- sort(unique(dat$Cluster))
	items_biopsy_sites <- sort(unique(dat$Biopsy_Site))

	ctr <- 1
	list.dm <- list()
	for(clust in items_clusters){
		for(biopsy_site in items_biopsy_sites){
			df <- data.frame(SampleID=dat$Sample_ID, Cluster=0, BiopsySite=0)
			df$Cluster[ which(dat$Cluster == clust) ] <- 1
			df$BiopsySite[ which(dat$Biopsy_Site == biopsy_site) ] <- 1

			# GET TABULAR MATRIX ---
			tbl <- as.matrix(table(df$Cluster, df$BiopsySite))
			tbl1 <- c(tbl[2,2], tbl[2,1], tbl[1,2], tbl[1,1])

			# CONFUSION MATRIX ---
			mat_conf <- matrix( tbl1, byrow=TRUE, nrow=2, ncol=2, dimnames=list(c("clust_yes","clust_no"),c("bsite_yes","bsite_no")))

			# FISHERS' EXACT TEST ---
			res_fishers_test <- fisher.test(x=mat_conf)

			list.dm[[ctr]] <- data.frame(Cluster=clust, 
											BiopsySite=biopsy_site, 
											Clust_yes_BSite_yes=tbl[2,2], 
											Clust_yes_BSite_no=tbl[2,1], 
											Clust_no_BSite_yes=tbl[1,2], 
											Clust_no_BSite_no=tbl[1,1],
											pvalue=res_fishers_test$p.value)
			ctr <- ctr + 1
		}
	}

	# AGGREGATE DATA ---
	dm <- do.call(rbind.data.frame, list.dm)
	rownames(dm) <- NULL

	return(dm)
}



########### FISHER'S EXACT TEST: ATAC-seq Clusters vs Biopsy Sites ###########
df <- subset(metadata, select=c("Sample_ID","Cluster","Biopsy_Site"))


### SUB-CLUSTERS 2A/2B ---
df2 <- subset(df, df$Cluster == "clust_2")
samplesids_2a <- c("DTB-080-PRO","DTB-009-BL","PR-095-BL","PR-056-BL","PR-120-BL","DTB-069-BL")
df2$Cluster <- ifelse(df2$Sample_ID %in% samplesids_2a, "sc_2A", "sc_2B")

#getFishersExactTest(dat=df2)

#> getFishersExactTest(dat=df2)
#  Cluster BiopsySite Clust_yes_BSite_yes Clust_yes_BSite_no Clust_no_BSite_yes  Clust_no_BSite_no      pvalue
#1   sc_2A       Bone                   1                  5                 16                  2 0.002704389
#2   sc_2A      Liver                   3                  3                  2                 16 0.078486731
#3   sc_2A Lymph Node                   2                  4                  0                 18 0.054347826
#4   sc_2B       Bone                  16                  2                  1                  5 0.002704389 #ASSOCIATION
#5   sc_2B      Liver                   2                 16                  3                  3 0.078486731
#6   sc_2B Lymph Node                   0                 18                  2                  4 0.054347826



### SUB-CLUSTERS 3A/3B ---
df3 <- subset(df, df$Cluster == "clust_3")
samplesids_3a <- c("DTB-085-BL","DTB-063-BL","DTB-080-BL","DTB-128-BL","DTB-167-PRO",
					"DTB-005-BL","DTB-259-BL","DTB-077-PRO","DTB-143-BL","DTB-019-PRO",
					"DTB-149-BL","DTB-119-PRO","DTB-035-BL","DTB-146-BL","DTB-101-BL",
					"DTB-167-BL","DTB-129-BL","DTB-034-BL","DTB-176-BL","DTB-055-PRO",
					"DTB-222-BL","DTB-124-BL","DTB-138-BL","PR-081-B")
df3$Cluster <- ifelse(df3$Sample_ID %in% samplesids_3a, "sc_3A", "sc_3B")

#getFishersExactTest(dat=df3)

#> getFishersExactTest(dat=df3)
#  Cluster BiopsySite Clust_yes_BSite_yes Clust_yes_BSite_no Clust_no_BSite_yes  Clust_no_BSite_no      pvalue
#1   sc_3A       Bone                  11                 12                  0                 15 0.00223446  #ASSOCIATION
#2   sc_3A      Liver                   0                 23                  3                 12 0.053935514
#3   sc_3A Lymph Node                   8                 15                 10                  5 0.095970498
#4   sc_3A      Other                   4                 19                  2                 13 1.000000000
#5   sc_3B       Bone                   0                 15                 11                 12 0.002234461
#6   sc_3B      Liver                   3                 12                  0                 23 0.053935514
#7   sc_3B Lymph Node                  10                  5                  8                 15 0.095970498
#8   sc_3B      Other                   2                 13                  4                 19 1.000000000

