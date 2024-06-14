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

dir.footprints <- file.path(dir.reproduce_data, "tf_footprints/mcrpc_subtypes")

### DEFINE FILES ---
file.rds_tf_hits <- file.path(dir.footprints, "data_footprints/tobias_mcrpc_subtypes_tf_hits.rds")

### SUBTYPES ---
subtypes <- c("ARpNEn","ARpNEp","ARlNEn","ARnNEp","ARnNEn")

### LOAD TF HITS---
list.tf_hits <- readRDS(file=file.rds_tf_hits)

### GET TF HITS SCORE ---
dat <- list.tf_hits$tf_footprints_hits_score

### ADD MOTIFID AND TF NAME ---
dat$motif_id <- unlist(lapply(stringr::str_split(dat$motif_prefix, "_"), function(x) x[1]))
dat$tf_name <- unlist(lapply(stringr::str_split(dat$motif_prefix, "_"), function(x) x[2]))

### RE-ORDER COLUMNS ----
dat <- subset(dat, select=c("motif_id","tf_name","ARpNEn","ARlNEn","ARpNEp","ARnNEp","ARnNEn"))


### WRITE OUTPUT ---
file.tbl <- file.path(dir.reproduce_tbl, "supplementary_table_s5.tsv")
write.table(dat, file.tbl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

cat("FILE GENERATED:", file.tbl, "\n", sep=" ")

