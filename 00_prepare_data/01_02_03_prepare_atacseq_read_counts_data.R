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
file.rds_atacseq_masterdata <- file.path(dir.reproduce_data, "wcdt_atacseq_masterdata.rds")
file.dat_mcrpc_raw <- file.path(dir.wrk, "analysis/02_explore_data/data/wcdt_atacseq_feature_count.tsv.gz")
file.dat_mcrpc_norm <- file.path(dir.wrk, "analysis/02_explore_data/data/wcdt_atacseq_feature_count_norm.tsv.gz")
file.annot_mcrpc <- file.path(dir.wrk, "analysis/02_explore_data/data/wcdt_atacseq_feature_annotation.tsv.gz")

file.dat_combined_raw <- file.path(dir.wrk, "analysis/03_pca_progression/data/feature_counts_combined_2021_10_31_exptWithCejas.tsv.gz")
file.dat_combined_norm <- file.path(dir.wrk, "analysis/03_pca_progression/data/feature_counts_combined_norm_2021_10_31_exptWithCejas.tsv.gz")
file.annot_combined <- file.path(dir.wrk, "analysis/03_pca_progression/data/feature_counts_combined_norm_2021_10_31_exptWithCejas_peak_annotation.tsv.gz")


### LOAD ATAC-SEQ MASTER DATA ---
list.rds_atacseq_masterdata <- readRDS(file=file.rds_atacseq_masterdata)



################### ATAC-SEQ READ COUNTS NORMALIZED: mCRPC --------
### LOAD FEATURE ANNOTATION DATA ---
annot_mcrpc <- data.table::fread(file=file.annot_mcrpc, sep="\t", header=TRUE, nThread=50, data.table=FALSE, verbose=FALSE)

### LOAD FEATURE COUNTS DATA: RAW ---
dat_mcrpc_raw <- data.table::fread(file=file.dat_mcrpc_raw, sep="\t", header=TRUE, nThread=50,  data.table=FALSE, verbose=FALSE)
dat_mcrpc_raw <- subset(dat_mcrpc_raw, select=c("feature_id", list.rds_atacseq_masterdata$metadata_mcrpc$Sample_Name))

### LOAD FEATURE COUNTS DATA: DESEQ2 NORMALIZED ---
dat_mcrpc_norm <- data.table::fread(file=file.dat_mcrpc_norm, sep="\t", header=TRUE, nThread=50,  data.table=FALSE, verbose=FALSE)
colnames(dat_mcrpc_norm)[1] <- "feature_id"
dat_mcrpc_norm <- subset(dat_mcrpc_norm, select=c("feature_id", list.rds_atacseq_masterdata$metadata_mcrpc$Sample_Name))

### ADD TO LIST ---
list.output_mcrpc <- list(feature_counts_raw=dat_mcrpc_raw,
                            feature_counts_norm=dat_mcrpc_norm,
                            feature_annotation=annot_mcrpc)

### SAVE OBJECT TO RDATA FILE ---
file.rds_atacseq_read_counts_mcrpc <- file.path(dir.reproduce_data, "atacseq_read_counts_mcrpc.rds")
saveRDS(object=list.output_mcrpc, file=file.rds_atacseq_read_counts_mcrpc)





################### ATAC-SEQ READ COUNTS NORMALIZED: COMBINED DATASET --------
### LOAD FEATURE ANNOTATION DATA ---
annot_combined <- data.table::fread(file=file.annot_combined, sep="\t", header=TRUE, nThread=50, data.table=FALSE, verbose=FALSE)

### LOAD FEATURE COUNTS DATA: RAW ---
dat_combined_raw <- data.table::fread(file=file.dat_combined_raw, sep="\t", header=TRUE, nThread=50,  data.table=FALSE, verbose=FALSE)
colnames(dat_combined_raw)[1] <- "feature_id"
dat_combined_raw <- subset(dat_combined_raw, select=c("feature_id", list.rds_atacseq_masterdata$metadata_combined$SAMPLE_ID))

### LOAD FEATURE COUNTS DATA: DESEQ2 NORMALIZED ---
dat_combined_norm <- data.table::fread(file=file.dat_combined_norm, sep="\t", header=TRUE, nThread=50,  data.table=FALSE, verbose=FALSE)
colnames(dat_combined_norm)[1] <- "feature_id"
dat_combined_norm <- subset(dat_combined_norm, select=c("feature_id", list.rds_atacseq_masterdata$metadata_combined$SAMPLE_ID))

### ADD TO LIST ---
list.output_combined <- list(feature_counts_raw=dat_combined_raw,
                            feature_counts_norm=dat_combined_norm,
                            feature_annotation=annot_combined)

### SAVE OBJECT TO RDATA FILE ---
file.rds_atacseq_read_counts_combined <- file.path(dir.reproduce_data, "atacseq_read_counts_combined.rds")
saveRDS(object=list.output_combined, file=file.rds_atacseq_read_counts_combined)
