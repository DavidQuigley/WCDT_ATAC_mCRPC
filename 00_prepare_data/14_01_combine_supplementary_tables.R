###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### DEFINE LIBRARIES ---
library("stringr")
library("data.table")
library("xlsx")

### DEFINE PATH ---
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/analysis/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

dir.footprints <- file.path(dir.reproduce_data, "tf_footprints/mcrpc_subtypes")

### DEFINE FILES ---
file.tbl1 <- file.path(dir.reproduce_tbl, "supplementary_table_s1.tsv")
file.tbl2 <- file.path(dir.reproduce_tbl, "supplementary_table_s2.tsv")
file.tbl3 <- file.path(dir.reproduce_tbl, "supplementary_table_s3.tsv")
file.tbl4 <- file.path(dir.reproduce_tbl, "supplementary_table_s4.tsv")
file.tbl5 <- file.path(dir.reproduce_tbl, "supplementary_table_s5.tsv")
file.tbl6 <- file.path(dir.reproduce_tbl, "supplementary_table_s6.tsv")
file.tbl_xls <- file.path(dir.reproduce_tbl, "supplementary_tables.xlsx")

### LOAD  DATA ---
dat1 <- data.table::fread(file.tbl1, data.table=FALSE, header=TRUE, stringsAsFactors=FALSE, nThread=1, showProgress=TRUE)
dat2 <- data.table::fread(file.tbl2, data.table=FALSE, header=TRUE, stringsAsFactors=FALSE, nThread=1, showProgress=TRUE)
dat3 <- data.table::fread(file.tbl3, data.table=FALSE, header=TRUE, stringsAsFactors=FALSE, nThread=1, showProgress=TRUE)
dat4 <- data.table::fread(file.tbl4, data.table=FALSE, header=TRUE, stringsAsFactors=FALSE, nThread=1, showProgress=TRUE)
dat5 <- data.table::fread(file.tbl5, data.table=FALSE, header=TRUE, stringsAsFactors=FALSE, nThread=1, showProgress=TRUE)
dat6 <- data.table::fread(file.tbl6, data.table=FALSE, header=TRUE, stringsAsFactors=FALSE, nThread=1, showProgress=TRUE)

### WRITE OUTPUT ---
xlsx::write.xlsx2(x=dat1, file=file.tbl_xls, sheetName="supplementary_table_s1", col.names=TRUE, row.names=FALSE, showNA=TRUE, append=FALSE)
xlsx::write.xlsx2(x=dat2, file=file.tbl_xls, sheetName="supplementary_table_s2", col.names=TRUE, row.names=FALSE, showNA=TRUE, append=TRUE)
xlsx::write.xlsx2(x=dat3, file=file.tbl_xls, sheetName="supplementary_table_s3", col.names=TRUE, row.names=FALSE, showNA=TRUE, append=TRUE)
xlsx::write.xlsx2(x=dat4, file=file.tbl_xls, sheetName="supplementary_table_s4", col.names=TRUE, row.names=FALSE, showNA=TRUE, append=TRUE)
xlsx::write.xlsx2(x=dat5, file=file.tbl_xls, sheetName="supplementary_table_s5", col.names=TRUE, row.names=FALSE, showNA=TRUE, append=TRUE)
xlsx::write.xlsx2(x=dat6, file=file.tbl_xls, sheetName="supplementary_table_s6", col.names=TRUE, row.names=FALSE, showNA=TRUE, append=TRUE)
