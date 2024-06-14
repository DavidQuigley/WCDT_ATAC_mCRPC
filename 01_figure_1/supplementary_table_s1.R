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
file.rds_rnaseq_ar_nepc_scores <- file.path(dir.reproduce_data, "wcdt_rnaseq_ar_nepc_scores.rds")

### LOAD ATAC-SEQ MASTER DATA ---
list.rds_atacseq_masterdata <- readRDS(file=file.rds_atacseq_masterdata)

### LOAD METADATA WGS ---
metadata_wgs <- subset(list.rds_atacseq_masterdata$stats_wgs_mcrpc, select=c("SampleID","Purity"))
colnames(metadata_wgs) <- c("Sample_ID","Tumor_Purity_WGS")


### LOAD RNA-SEQ AR/NEPC SCORE DATA ---
dat_ar_nepc_scores <- readRDS(file=file.rds_rnaseq_ar_nepc_scores)
colnames(dat_ar_nepc_scores) <- c("Sample_ID","ARscore","NEscore")

### GET METADATA ---
metadata <- list.rds_atacseq_masterdata$metadata_mcrpc

### BIOPSY SITE ---
metadata$Biopsy_Site <- metadata$Biopsy_Site_detail

### SUBTYPE ---
metadata$mCRPC_Subtype <- stringr::str_replace_all(stringr::str_replace_all(metadata$mCRPC_Subtype, "_", ""), "ARLNE-", "ARlowNE-")

### EXTRACT COLUMNS ---
items <- c("Sample_ID","Sample_Name","Patient_ID","Processed_Core","Biopsy_Site","mCRPC_Subtype","ATACseq_accession","RNAseq_accession",
                    "ATACseq_Mapped_Reads","ATACseq_Peak_Counts","ATACseq_FRiP","ATACseq_TSS_Enrich_Score")

dat <- subset(metadata, select=items)

### MERGE DATA ---
dat <- merge(dat, metadata_wgs, by="Sample_ID", all=TRUE)
dat <- merge(dat, dat_ar_nepc_scores, by="Sample_ID", all=TRUE)

### WRITE OUTPUT ---
file.tbl <- file.path(dir.reproduce_tbl, "supplementary_table_s1.tsv")
write.table(dat, file.tbl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

cat("FILE GENERATED:", file.tbl, "\n", sep=" ")
