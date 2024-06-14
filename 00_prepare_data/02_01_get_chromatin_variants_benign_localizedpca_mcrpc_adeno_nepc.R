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

dir.diffpeaks <- file.path(dir.wrk, "analysis/03_pca_progression/data")

### DEFINE FILES ---
file.diffpeaks_1 <- file.path(dir.diffpeaks, "diffpeaks_benign_vs_pca.tsv.gz")
file.diffpeaks_2 <- file.path(dir.diffpeaks, "diffpeaks_benign_vs_mcrpc_adeno.tsv.gz")
file.diffpeaks_3 <- file.path(dir.diffpeaks, "diffpeaks_benign_vs_mcrpc_nepc.tsv.gz")
file.diffpeaks_4 <- file.path(dir.diffpeaks, "diffpeaks_pca_vs_mcrpc_adeno.tsv.gz")
file.diffpeaks_5 <- file.path(dir.diffpeaks, "diffpeaks_pca_vs_mcrpc_nepc.tsv.gz")
file.diffpeaks_6 <- file.path(dir.diffpeaks, "diffpeaks_mcrpc_adeno_vs_nepc.tsv.gz")

file.clust_promoter <- file.path(dir.diffpeaks, "combinedata_heatmap_atacseq_feature_cluster_promoter.tsv")
file.clust_intron <- file.path(dir.diffpeaks, "combinedata_heatmap_atacseq_feature_cluster_intron.tsv")
file.clust_intergenic <- file.path(dir.diffpeaks, "combinedata_heatmap_atacseq_feature_cluster_intergenic.tsv")

### COMPILE FILE PATHS ---
files.diffpeaks <- c(file.diffpeaks_1, file.diffpeaks_2, file.diffpeaks_3, file.diffpeaks_4, file.diffpeaks_5, file.diffpeaks_6)
y_analysis <- c("benign_vs_pca", "benign_vs_mcrpc_adeno", "benign_vs_mcrpc_nepc", "pca_vs_mcrpc_adeno", "pca_vs_mcrpc_nepc", "mcrpc_adeno_vs_nepc")
names(files.diffpeaks) <- y_analysis

### LOAD DATA ---
list.diffpeaks <- list()
for(analysis in y_analysis){
    list.diffpeaks[[analysis]] <- data.table::fread(file=files.diffpeaks[analysis], sep="\t", header=TRUE, nThread=50, data.table=FALSE, verbose=FALSE)
}


### FUNCTION: filterData() ---
filterData <- function(dat, pvalue.threshold=0.01, logfc.threshold=1){
    # DIFFERENTIAL ACCESSIBILITY ---
    d <- subset(dat, dat$padj <= pvalue.threshold)
    d <- subset(d, (d$log2FoldChange >= logfc.threshold) | (d$log2FoldChange <= -logfc.threshold) )
    d$fc_direction <- ifelse( sign(d$log2FoldChange) == 1, "UP", "DOWN")

    return(d)
}


### FILTER DATA FOR EACH COMPARISON ---
list.diffpeaks_pass <- lapply(list.diffpeaks, function(x) filterData(dat=x, pvalue.threshold=0.01, logfc.threshold=1) )


### ADD UP/DN TAG ----
for(analysis in y_analysis){
    df <- list.diffpeaks_pass[[analysis]]

    # TAG PHENOTYPE BY UP OR DOWN ACCESSIBILITY ---    
    df$UP <- apply(df, 1, function(x) ifelse(x[17] == "UP", stringr::str_split(analysis, "_vs_")[[1]][2], stringr::str_split(analysis, "_vs_")[[1]][1]) )
    df$DN <- apply(df, 1, function(x) ifelse(x[17] != "UP", stringr::str_split(analysis, "_vs_")[[1]][2], stringr::str_split(analysis, "_vs_")[[1]][1]) )

    # REPLACE ---
    df$UP <- stringr::str_replace_all(df$UP, "nepc", "mcrpc_nepc")
    df$DN <- stringr::str_replace_all(df$DN, "nepc", "mcrpc_nepc")

    # STORE DATA ---
    list.diffpeaks_pass[[analysis]] <- df
}


#################################################################################################################################################################
### AGGREGATE DATA ---
df <- do.call(rbind.data.frame,  lapply(list.diffpeaks_pass, function(x) subset(x, select=c("FeatureID","Feature","Gene","GeneType","distanceToTSS","UP","DN") ) ) )
rownames(df) <- NULL

### GET UNIQUE SET OF FETUREIDS ---
dm <- subset(df, select=c("FeatureID","Feature","Gene","GeneType","distanceToTSS"))
dm <- dm[!duplicated(dm),]

### ADD COLUMNS ---
dm$benign <- 0
dm$pca <- 0     # localized PCa
dm$mcrpc_adeno <- 0
dm$mcrpc_nepc <- 0

### TAG FEATURES WITH DIFF PEAKS UP ---
phenotypes <- c("benign","pca","mcrpc_adeno","mcrpc_nepc")

for(phenotype in phenotypes){
    feature_ids <- unique(df$FeatureID[which(df$UP == phenotype)])
    dm[which(dm$FeatureID %in% feature_ids), phenotype] <- 1
}


#################################################################################################################################################################
### BASED ON HEATMAP CLUSTERING OF ATAC-SEQ READCOUNTS -----
### LOAD CLUSTER DATA ---
dat_clust_promoter <- data.table::fread(file=file.clust_promoter, sep="\t", header=TRUE, nThread=1, data.table=FALSE, verbose=FALSE)
dat_clust_intron <- data.table::fread(file=file.clust_intron, sep="\t", header=TRUE, nThread=1, data.table=FALSE, verbose=FALSE)
dat_clust_intergenic <- data.table::fread(file=file.clust_intergenic, sep="\t", header=TRUE, nThread=1, data.table=FALSE, verbose=FALSE)

### MERGE DATA ---
dat_clust <- rbind(dat_clust_promoter, dat_clust_intron, dat_clust_intergenic)



#################################################################################################################################################################
### ADD TO LIST ---
list.output <- list(atacseq_diffpeaks=list.diffpeaks,
                    atacseq_diffpeaks_pass=list.diffpeaks_pass,
                    atacseq_diffpeaks_up_open=dm,
                    atacseq_diffpeaks_cluster_features=dat_clust)


### SAVE OBJECT TO RDATA FILE ---
file.rds_atacseq_chromatin_variants <- file.path(dir.reproduce_data, "atacseq_chromatin_variants_benign_localizedpca_mcrpc_adeno_nepc.rds")
saveRDS(object=list.output, file=file.rds_atacseq_chromatin_variants)
