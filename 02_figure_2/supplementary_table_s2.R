###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### DEFINE PATH ---
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

### DEFINE FILES ----
file.rds_genelist <- file.path(dir.reproduce_data, "genelist.rds")



### LOAD GENELIST DATA ---
list.rds_genelist <- readRDS(file=file.rds_genelist)
genelist.ar <- list.rds_genelist$msigdb_pathways$hallmark$HALLMARK_ANDROGEN_RESPONSE 
genes_nepc <- list.rds_genelist$genes_nepc


### GET MAX LENGTH OF GENES ---
n_maxlength <- max( length(genelist.ar), length(genes_nepc) )

### GET VECTORS OF EQAL LENGTH ---
vec_ar <- c(genelist.ar, rep("", n_maxlength - length(genelist.ar)))  
vec_ne <- c(genes_nepc, rep("", n_maxlength - length(genes_nepc)))  

### COMBINE TWO VECTORS ---
df <- cbind(vec_ar, vec_ne)
colnames(df) <- c("genelist_AR_score","genelist_NE_score")

### WRITE OUTPUT ---
file.tbl <- file.path(dir.reproduce_tbl, "supplementary_table_s2.tsv")
write.table(df, file.tbl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
