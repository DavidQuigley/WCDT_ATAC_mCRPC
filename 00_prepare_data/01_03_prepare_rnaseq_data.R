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
file.rnaseq_wcdt <- file.path("/data1/datasets_1/human_prostate_WCDT/rna_lasercapture_2018/processed/counts_genocode28/Counts_TPM_Genes.txt")
file.rnaseq_promote <- file.path("/data1/datasets_1/human_prostate_WCDT/rna_promote/processed/counts_gencode28/Counts_TPM_Genes.txt")
file.annotation_gtf <- file.path(dir.wrk, "reference/annotation_data/gene_annotation/gencode_v28/gencode.v28.annotation.gtf.gz")


### FUNCTION: parseGTF() --
parseGTF <- function(file.gtf, feature.type){
    # CHROMOSOMES ---
    chromosome <- paste("chr", c(1:22,"X","Y"), sep="")

    # LOAD GTF FILE ---
    annot <- data.table::fread(file=file.gtf, sep="\t", header=FALSE, nThread=50, data.table=TRUE, stringsAsFactors=FALSE, skip=5, verbose=FALSE)
    annot <- subset(annot, annot$V1 %in% chromosome)

    # PARSE GTF COLUMNS ---
    annot <- subset(annot, annot$V3 == "gene")
    list.annot <- noquote(stringr::str_split(annot$V9, "; "))
    annot$EnsemblID <- unlist(lapply(stringr::str_split(unlist(lapply(list.annot, function(x) gsub('"', '', x[stringr::str_detect(x, "gene_id")]))), " "), function(x) x[2]))
    annot$Gene <- unlist(lapply(stringr::str_split(unlist(lapply(list.annot, function(x) gsub('"', '', x[stringr::str_detect(x, "gene_name")]))), " "), function(x) x[2]))
    annot$GeneType <- unlist(lapply(stringr::str_split(unlist(lapply(list.annot, function(x) gsub('"', '', x[stringr::str_detect(x, "gene_type")]))), " "), function(x) x[2]))

    # TRIM DATA ---
    #annot <- subset(annot, annot$GeneType == "protein_coding")   
    items <- c("V1","V4","V5","Gene","EnsemblID","V7","GeneType")
    df <- subset(annot, select=items)
    colnames(df) <- c("Chrom","Start","End","Gene","EnsemblID","Strand","GeneType")
    
    return(df)
}


### LOAD ATAC-SEQ MASTER DATA ---
list.rds_atacseq_masterdata <- readRDS(file=file.rds_atacseq_masterdata)
metadata <-list.rds_atacseq_masterdata$metadata_mcrpc

### EXTRACT GENE ANNOTATION ---
annot <- parseGTF(file.gtf=file.annotation_gtf, feature.type="gene")


### LOAD RNA-SEQ GENE EXPRESSION DATASET ---
dat_wcdt <- data.table::fread(file=file.rnaseq_wcdt, sep="\t", header=TRUE, nThread=50,  data.table=TRUE, verbose=FALSE)
dat_promote <- data.table::fread(file=file.rnaseq_promote, sep="\t", header=TRUE, nThread=50,  data.table=TRUE, verbose=FALSE)
dat_wcdt$gene_name <- NULL
dat_promote$gene_name <- NULL


### MERGE EXPRESSION DATA BY FEATURE IDS ---
expr <- merge(dat_wcdt, dat_promote, by="FEATURE_ID")

### EXTRACT EXPRESSION DATA FOR ATAC-SEQ SAMPLES ---
expr_atac <- subset(expr, select=c("FEATURE_ID", metadata$SampleID_RNAseq))




#### GET COMMON FEATURE IDS ---
featureids <- intersect(expr_atac$FEATURE_ID, annot$EnsemblID)

### EXTRACT EXPRESSION BY FEATUREIDS ---
expr_atac <- subset(expr_atac, expr_atac$FEATURE_ID %in% featureids)
annot <- subset(annot, annot$EnsemblID %in% featureids)

### ORDER EXPRESSION BY FEATUTEIDS ---
expr_atac <- expr_atac[match(annot$EnsemblID, expr_atac$FEATURE_ID),]

### ADD TO LIST ---
list.output <- list(gene_annotation=annot, gene_expression=expr_atac)


### SAVE OBJECT TO RDATA FILE ---
file.rds_rnaseq <- file.path(dir.reproduce_data, "wcdt_rnaseq_atacseq_gene_expression.rds")
saveRDS(object=list.output, file=file.rds_rnaseq)
