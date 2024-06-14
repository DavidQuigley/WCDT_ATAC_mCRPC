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

dir.input <- file.path(dir.reproduce_data, "tfpeaks2gene/input")

### DEFINE FILES ---
file.rds_atacseq_masterdata <- file.path(dir.reproduce_data, "wcdt_atacseq_masterdata.rds")
file.rds_rnaseq <- file.path(dir.reproduce_data, "wcdt_rnaseq_atacseq_gene_expression.rds")
file.rds_atacseq_read_counts_mcrpc <- file.path(dir.reproduce_data, "atacseq_read_counts_mcrpc.rds")
file.rds_peaks_footprints <- file.path(dir.reproduce_data, "wcdt_atacseq_peak_footprint_overlap.rds")

### CHROMOSOME ---
chromosomes <- paste("chr", c(1:22, "X","Y"), sep="")


##################################################################################################################
### LOAD ATAC-SEQ MASTER DATA ---
list.rds_atacseq_masterdata <- readRDS(file=file.rds_atacseq_masterdata)
metadata <- list.rds_atacseq_masterdata$metadata_mcrpc




##################################################################################################################
########## PREPARE ATAC-SEQ DATA ---------------------------------------------------------------------------------

### LOAD ATAC-SEQ:TF-FOOTPRINTS DATA ---
list.rds_peaks_footprints <- readRDS(file=file.rds_peaks_footprints)
motifs <- names(list.rds_peaks_footprints)

### MOTIF FOOTPRINTS DATA TO PEAKS ---
list.annot_peaks <- list()
for(motif_prefix in motifs){
    df <- do.call(rbind.data.frame, list.rds_peaks_footprints[[motif_prefix]])
    rownames(df) <- NULL
    list.annot_peaks[[motif_prefix]] <- df
}

annot_peaks <- do.call(rbind.data.frame, list.annot_peaks)
rownames(annot_peaks) <- NULL
annot_peaks$chr <- unlist(lapply(stringr::str_split(annot_peaks$peak_position, "_"), function(x) x[1]))


### LOAD ATAC-SEQ DATA ---
list.rds_atacseq_read_counts_mcrpc <- readRDS(file=file.rds_atacseq_read_counts_mcrpc)
dat_atacseq <- list.rds_atacseq_read_counts_mcrpc$feature_counts_norm
rownames(dat_atacseq) <- dat_atacseq$feature_id
dat_atacseq$feature_id <- NULL
dat_atacseq <- subset(dat_atacseq, select=metadata$Sample_Name)
colnames(dat_atacseq) <- metadata$Sample_ID



### EXTRACT DATA BY CHROMOSOME ---
for(chr in chromosomes){
    cat("ATAC-SEQ START:", chr, "\n", sep=" ")

    annot_peaks_chr <- annot_peaks[which(annot_peaks$chr == chr),]
    dat_atacseq_chr <- dat_atacseq[which(rownames(dat_atacseq) %in% unique(annot_peaks_chr$peak_position)),]
    mat_atacseq <- as.matrix(t(dat_atacseq_chr))

    # WRITE OUTPUT ----
    file.atacseq_out <- file.path(dir.input, paste("atacseq_peaks2gene_input_", chr,  ".tsv", sep=""))
    write.table(mat_atacseq, file.atacseq_out, sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)

    # COMPRESS ---
    cmd <- paste("gzip", file.atacseq_out, sep=" ")
    system(cmd)

    cat("ATAC-SEQ PROCESSED:", chr, ncol(mat_atacseq), "\n", sep=" ")
}





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


### EXTRACT RNA-SEQ DATA BY CHROMOSOME ---
for(chr in chromosomes){
    cat("RNA-SEQ START:", chr, "\n", sep=" ")

    genes <- annot_gene$Gene[which(annot_gene$Chrom == chr)]
    expr_chr <- expr[which(rownames(expr) %in% genes),]
    mat_rnaseq <- as.matrix(t(expr_chr))

    # WRITE OUTPUT ----
    file.rnaseq_out <- file.path(dir.input, paste("rnaseq_peaks2gene_input_", chr, ".tsv", sep=""))
    write.table(mat_rnaseq, file.rnaseq_out, sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)

    # COMPRESS ---
    cmd <- paste("gzip", file.rnaseq_out, sep=" ")
    system(cmd)

    cat("RNA-SEQ PROCESSED:", chr, ncol(mat_rnaseq), "\n", sep=" ")
}
