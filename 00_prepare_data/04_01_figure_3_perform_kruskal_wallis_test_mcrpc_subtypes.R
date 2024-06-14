###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### DEFINE LIBRARIES ---
library("stringr")
library("data.table")
library("dplyr")
library("reshape2")
library("foreach")
library("doParallel")

### DEFINE PATH ---
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/analysis/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

### DEFINE FILES ----
file.rds_atacseq_masterdata <- file.path(dir.reproduce_data, "wcdt_atacseq_masterdata.rds")
file.rds_atacseq_read_counts_mcrpc <- file.path(dir.reproduce_data, "atacseq_read_counts_mcrpc.rds")



### FUNCTION: prepareData() ---
prepareData <- function(mat, metadata){
    # RESHAPE DATA ---
    df <- reshape2::melt(t(as.matrix(mat)))
    colnames(df) <- c("SampleID","FeatureID","Value")

    # TO CHARACTER ---
    df$SampleID <- as.character(df$SampleID)
    df$FeatureID <- as.character(df$FeatureID)

    # GET UNIQUE SUBTYPES ---
    subtypes <- c("AR+_NE-","ARL_NE-","AR+_NE+","AR-_NE+","AR-_NE-")

    # ADD GROUPS ---
    df$Group <- NA
    for(subtype in subtypes){
        sampleids_subtype <- metadata$Sample_ID[which(metadata$mCRPC_Subtype == subtype)]
        df$Group[which(df$SampleID %in% sampleids_subtype)] <- subtype
    }

    return(df)
}


### FUNCTION: calcKruskalWallisTest() --------
calcKruskalWallisTest <- function(df, feature_name, n_cores=50){
	# GET ATTRIBUTES ---
    colnames(df)[which(colnames(df) == feature_name)] <- "Feature"
	features <- as.character(unique(df$Feature))
	nfeatures <- length(features)

	# Declate Cluster ----
	cl <- makeCluster(n_cores)
	registerDoParallel(cl)

	# LOOP FOR EACH FEATURE ---
	dm <- foreach(i = 1:nfeatures, .inorder=FALSE, .combine=rbind, .errorhandling="remove", .verbose=TRUE) %dopar%{	
			    d <- subset(df, df$Feature == features[i])
				kwtest <- kruskal.test(Value ~ Group, data = d, na.action=na.omit)

				dm <- data.frame(Feature=features[i],
								statistic=as.numeric(kwtest$statistic),
								pvalue=kwtest$p.value)
			}

	stopCluster(cl)

	dm <- dm[order(dm$pvalue, decreasing=FALSE),]
	dm$fdr <- p.adjust(dm$pvalue, method = "BH", n = length(dm$pvalue))
	
	return(dm)
}




######################################################################
### LOAD ATAC-SEQ MASTER DATA ---
list.rds_atacseq_masterdata <- readRDS(file=file.rds_atacseq_masterdata)
metadata <- list.rds_atacseq_masterdata$metadata_mcrpc



### LOAD ATAC-SEQ READ COUNTS DATA ---
list.rds_atacseq_read_counts_mcrpc <- readRDS(file=file.rds_atacseq_read_counts_mcrpc)


### GET NORMALIZED ATAC-SEQ READCOUNTS DATA ---
mat <- list.rds_atacseq_read_counts_mcrpc$feature_counts_norm
rownames(mat) <- mat$feature_id
mat$feature_id <- NULL
colnames(mat) <- metadata$Sample_ID




######################################################################
### PREPARE DATA ---
df <- prepareData(mat, metadata, n_cores=50)


### KRUSKAL WALLIS TEST ---
dm_kwtest <- calcKruskalWallisTest(df, feature_name="FeatureID")


### WRITE OUTPUT ---
file.output <- file.path(dir.reproduce_data, "wcdt_atacseq_mcrpc_subtypes_KruskalWallisTest_summary.tsv")
write.table(dm_kwtest, file.output, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

### COMPRESS ---
cmd <- paste("gzip", file.output, sep=" ")
system(cmd)
