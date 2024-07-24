###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### LOAD LIBRARIES ---
library("stringr")
library("data.table")
library("tibble")
library("dplyr")


### DEFINE PATH ---
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

dir.tfpeaks2gene <- file.path(dir.reproduce_data, "tfpeaks2gene")
dir.footprints <- file.path(dir.reproduce_data, "tf_footprints/mcrpc_subtypes/data_footprints") 

### DEFINE FILES ---
file.tf2genes <- file.path(dir.tfpeaks2gene, "wcdt_tf2genes_r_0p4_pval_0p05_distTSS_500kb.tsv.gz")
file.tf_hits_rds <- file.path(dir.footprints, "tobias_mcrpc_subtypes_tf_hits.rds")

### SUBTYPES --
subtypes <- c("ARpNEn","ARlNEn","ARpNEp","ARnNEp","ARnNEn")
subtypes_final <- c("AR+NE-","AR+NE+","ARlowNE-","AR-NE+","AR-NE-")


###########################################################################################################################
### FUNCTION: get_tf2gene_pairs() ---
get_tf2gene_pairs <- function(dat_tf2genes, dat_tfsubtypes, subtype){
    # EXTRACT TF BY SUBTYPE ---
    index_subtype <- which( colnames(dat_tfsubtypes) == subtype )
    index_tf <- which( dat_tfsubtypes[,index_subtype] >= 2 ) 
    dat_tfsubtypes <- dat_tfsubtypes[index_tf,]
    tfs_subtype <- dat_tfsubtypes$motif_prefix

    # EXTRACT TF2GENES BY TF ---
    dat_tf2genes_sub <- subset(dat_tf2genes, dat_tf2genes$motif_prefix %in% tfs_subtype)

    # EXTRACT TF-GENE PAIRS ---
    #df <- subset(dat_tf2genes_sub, select=c("motif_prefix","Gene"))
    #df <- df[!duplicated(df),]

    return(dat_tf2genes_sub)
}

### FUNCTION: prepareData() ---
prepareData <- function(df){
    # GET COUNT TABLE ---
    dm <- tibble::tibble(df) %>% 
            dplyr::select(motif_prefix, Gene) %>%
            dplyr::distinct(motif_prefix, Gene, .keep_all = TRUE) %>%
            dplyr::count(motif_prefix, sort=TRUE, name="Freq")

    dm$tf_name <- unlist(lapply(str_split(dm$motif_prefix, "_"), function(x) paste(x[2:length(x)], collapse="_") ))
    dm$Rank <- c(1:nrow(dm))    

    return(dm)
}

###########################################################################################################################
### LOAD TF2GENES DATA ---
dat_tf2genes <- data.table::fread(file=file.tf2genes, sep="\t", header=TRUE, nThread=50, data.table=FALSE, verbose=FALSE)

### LOAD TF HITS DATA ---
list.tf_hits <- readRDS(file=file.tf_hits_rds)
dat_tfsubtypes <- list.tf_hits$tf_footprints_hits_refined


###########################################################################################################################
### GET TF2GENE PAIRS PER SUBTYPE ----
list.tf2gene <- list()
list.tf2gene_stats <- list()
for(subtype in subtypes){
    list.tf2gene[[subtype]] <- get_tf2gene_pairs(dat_tf2genes, dat_tfsubtypes, subtype)
    list.tf2gene_stats[[subtype]] <- prepareData(df=list.tf2gene[[subtype]])

    cat("PROCESSED:", subtype, "\n", sep=" ")
}


### PREPARE DATA FOR TABLE ---
list.tbl <- list()
for(subtype in subtypes){
    dat <- list.tf2gene_stats[[subtype]]

    list.tbl[[subtype]] <- data.frame(subtype=subtype,
                                        motif_id=unlist(lapply(stringr::str_split(dat$motif_prefix, "_"), function(x) x[1])),
                                        tf_name=unlist(lapply(stringr::str_split(dat$motif_prefix, "_"), function(x) x[2])),
                                        number_of_target_genes=dat$Freq)   
}

### AGGREGATE DATA ---
df <- do.call(rbind.data.frame, list.tbl)
rownames(df) <- NULL

### CHANGE mCRPC SUBTYPE ANNOTATION ---
for(i in 1:length(subtypes)){
    df$subtype[which(df$subtype == subtypes[i])] <- subtypes_final[i]
}



### WRITE OUTPUT ---
file.tbl <- file.path(dir.reproduce_tbl, "supplementary_table_s6.tsv")
write.table(df, file.tbl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

cat("FILE GENERATED:", file.tbl, "\n", sep=" ")

