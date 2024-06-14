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

dir.permute <- file.path(dir.reproduce_data, "tfpeaks2gene/permutation")

### DEFINE FILES ---
file.rds_atacseq_masterdata <- file.path(dir.reproduce_data, "wcdt_atacseq_masterdata.rds")
file.rds_rnaseq <- file.path(dir.reproduce_data, "wcdt_rnaseq_atacseq_gene_expression.rds")

### CHROMOSOME ---
chromosomes <- paste("chr", c(1:22, "X","Y"), sep="")


##################################################################################################################
### LOAD ATAC-SEQ MASTER DATA ---
list.rds_atacseq_masterdata <- readRDS(file=file.rds_atacseq_masterdata)
metadata <- list.rds_atacseq_masterdata$metadata_mcrpc



########## PREPARE RNA-SEQ EXPRESSION DATA -----------------------------------------------------------------------------------
### LOAD RNA-SEQ GENE EXPRESSION OBJECT: TPM ---
cat(format(Sys.time(), "%b %d %X"), "LOADAING R DATA-OBJECT: RNA-seq EXPRESSION ... ", "\n", sep=" ")
    list.rnaseq <- readRDS(file=file.rds_rnaseq)

    # ANNOTATION ---
    gannot <- list.rnaseq$gene_annotation
    gannot <- subset(gannot, gannot$GeneType == "protein_coding")
    gannot <- gannot[-which(stringr::str_detect(gannot$EnsemblID, "PAR_Y")),]

    # EXPRESSION ---
    expr <- as.data.frame(list.rnaseq$gene_expression)
    expr <- subset(expr, expr$FEATURE_ID %in% gannot$EnsemblID)
    rownames(expr) <- gannot$Gene
    expr$FEATURE_ID <- NULL
    colnames(expr) <- stringr::str_replace(colnames(expr), "RP", "BL")

cat(format(Sys.time(), "%b %d %X"), "DONE! ", "\n", sep=" ")

### RE-ARRANGE EXPRESSION DATA ---
expr <- expr[order(rownames(expr), decreasing=FALSE),]
expr <- subset(expr, select=metadata$Sample_ID)

### LOG2 TRANSFORM [log2(TPM + 1)]---
expr <- log2(expr + 1)

########## -------------------------------------------------------------------------------------------------------------------

### REMOVE GENES WITH STANDARD DEVIATION = 0 ---
std_dev <- apply(expr, 1, sd)
del.index_1 <- which(std_dev == 0)
expr <- expr[-del.index_1,]

### REMOVE GENES WITH VALUE TPM <= 0 IN > 25% OF SAMPLES ---
n <- round(ncol(expr) * (25/100), 0)
ln_0 <- apply(expr, 1, function(x) length(which(x <= log2(0+1) )))
del.index_2 <- which(ln_0 >= n)
expr <- expr[-del.index_2,]

### GET GENE ANNOTATION ---
annot_gene <- subset(gannot, gannot$Gene %in% rownames(expr))




### GET RANDOMIZED GENE ORDER ---
genelist_original <- rownames(expr)
n <- 1000
list.genes <- list()
set.seed(12345)
for(i in 1:n){
    list.genes[[i]] <- sample(genelist_original, length(genelist_original))
}

### EXTRACT INDICES OF RANDOMIZE GENE LABELS ---
list.permutated_gene_index <- list()
for(chr in chromosomes){
    cat("START:", chr, "\n", sep=" ")
    genes_chr <- annot_gene$Gene[which(annot_gene$Chrom == chr)]

    list.gene_index <- list()
    for(j in 1:n){
        list.gene_index[[j]] <- which(list.genes[[j]] %in% genes_chr)
    }

    list.permutated_gene_index[[chr]] <- list.gene_index
    
    cat("PROCESSED:", chr, "\n", sep=" ")
}


### PREPARE DATA INTO LIST ---
list.perm_expr <- list(expression=t(expr), permutated_genes=list.genes, permutated_gene_index=list.permutated_gene_index)

### SAVE OBJECT TO RDATA FILE ---
file.perm_rds <- file.path(dir.permute, "wcdt_rnaseq_permutation.rds")
saveRDS(object=list.perm_expr, file=file.perm_rds)
