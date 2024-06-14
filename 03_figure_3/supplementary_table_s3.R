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

### DEFINE FILES ----
file.rds_atacseq_read_counts_mcrpc <- file.path(dir.reproduce_data, "atacseq_read_counts_mcrpc.rds")
file.kwtest <- file.path(dir.reproduce_data, "wcdt_atacseq_mcrpc_subtypes_KruskalWallisTest_summary.tsv.gz")


### LOAD ATAC-SEQ READ COUNTS DATA ---
list.rds_atacseq_read_counts_mcrpc <- readRDS(file=file.rds_atacseq_read_counts_mcrpc)

### LOAD ATAC-SEQ PEAK ANNOTATION DATA ---
annot <- list.rds_atacseq_read_counts_mcrpc$feature_annotation

### LOAD KRUSKAL-WALLIS TEST DATA ---
dat_kwtest <- data.table::fread(file=file.kwtest, sep="\t", header=TRUE, nThread=50, data.table=FALSE, verbose=FALSE)
colnames(dat_kwtest)[1] <- "FeatureID"
dat_kwtest <- merge(dat_kwtest, annot, by="FeatureID")
dat_kwtest <- dat_kwtest[order(dat_kwtest$pvalue, decreasing=FALSE),]
dat_kwtest <- subset(dat_kwtest, dat_kwtest$pvalue <= 0.001)
rownames(dat_kwtest) <- NULL

dat_kwtest$seqnames <- NULL
dat_kwtest$start <- NULL
dat_kwtest$end <- NULL

### COLNAMES ---
colnames(dat_kwtest) <- c("feature_id","statistic","pvalue","fdr","gene","gene_id","transcript_id","gene_type","distance_to_tSS","feature")

### WRITE OUTPUT ---
file.tbl <- file.path(dir.reproduce_tbl, "supplementary_table_s3.tsv")
write.table(dat_kwtest, file.tbl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

cat("FILE GENERATED:", file.tbl, "\n", sep=" ")

