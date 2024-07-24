###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### DEFINE LIBRARIES ---
library("stringr")

### DEFINE PATH ---
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

dir.tfpeaks2gene <- file.path(dir.reproduce_data, "tfpeaks2gene")
dir.corr <- file.path(dir.reproduce_data, "tfpeaks2gene/correlation") 


dir.scripts <- file.path(dir.reproduce, "scripts_temp")  # REPLACE PATH
dir.logs <- tempdir() # REPLACE PATH

### DEFINE FILES ----
file.rds_atacseq_read_counts_mcrpc <- file.path(dir.reproduce_data, "atacseq_read_counts_mcrpc.rds")

### CHROMOSOME ---
chromosomes <- paste("chr", c(1:22, "X","Y"), sep="")

### PARAMS ---
list.params <- list()
list.params[[1]] <- list(r=0.4, pval=0.05, correct_pval=FALSE, key="r_0p4_pval_0p05")
list.params[[2]] <- list(r=0.5, pval=0.05, correct_pval=FALSE, key="r_0p5_pval_0p05")
list.params[[3]] <- list(r=0.5, pval=0.01, correct_pval=FALSE, key="r_0p5_pval_0p01")
list.params[[4]] <- list(r=0.5, pval=0.05, correct_pval=TRUE, key="r_0p5_fdr_0p05")
list.params[[5]] <- list(r=0.5, pval=0.01, correct_pval=TRUE, key="r_0p5_fdr_0p01")
list.params[[6]] <- list(r=0.6, pval=0.05, correct_pval=FALSE, key="r_0p6_pval_0p05")
list.params[[7]] <- list(r=0.7, pval=0.05, correct_pval=FALSE, key="r_0p7_pval_0p05")
list.params[[8]] <- list(r=0.6, pval=0.05, correct_pval=TRUE, key="r_0p6_fdr_0p05")
list.params[[9]] <- list(r=0.7, pval=0.05, correct_pval=TRUE, key="r_0p7_fdr_0p05")



### FUNCTION: load_data() ---
load_data <- function(file.dat){
    # LOAD DATA ---
    dat <- data.table::fread(file=file.dat, sep="\t", header=TRUE, nThread=60, data.table=FALSE, verbose=FALSE)

    # FILL VALUE FOR PVALUE == 0
    index_0 <- which(dat$pvalue == 0)

    if(length(index_0) != 0){
        dat$pvalue[which(dat$pvalue == 0)] <- min(setdiff(dat$pvalue, 0))
    }

    # GET FDR PVALUE ---
    dat$fdr <- p.adjust(p=dat$pvalue, method="fdr", n=length(dat$pvalue))

    # LOG TRANSFORM PVALUES
    dat$nlogp <- -log10(dat$pvalue)
    dat$nlogq <- -log10(dat$fdr)

    return(dat)
}

### FUNCTION: filter_data() ---
filter_data <- function(list.dat, r, pval, correct_pval=TRUE){
    if(correct_pval){
        list.d1 <- lapply(list.dat, function(x) x[which(x$fdr <= pval),])
        list.d2 <- lapply(list.d1, function(x) x[which(abs(x$Corr) >= r),])
    }else{
        list.d1 <- lapply(list.dat, function(x) x[which(x$pvalue <= pval),])
        list.d2 <- lapply(list.d1, function(x) x[which(abs(x$Corr) >= r),])
    }

    # AGGREAGATE DATA ---
    df <- do.call(rbind.data.frame, list.d2)
    rownames(df) <- NULL

    return(df)
}


### GET PEAKS2GENE CORRELATION BY CHROMOSOME ---
list.dat <- list()
for(i in 1:length(chromosomes)){
    chr <- chromosomes[i]

    cat(format(Sys.time(), "%X"), "START:", chr, "\n", sep=" ")  

    # DEFINE FILE ---
    file.dat <- file.path(dir.corr, paste("pvalue_peaks2gene_corr_", chr,  ".tsv.gz", sep=""))

    # STORE DATA ---
    list.dat[[chr]] <- load_data(file.dat)

    cat(format(Sys.time(), "%X"), "PROCESSED:", chr, "\n", sep=" ")  
}

### FILTER DATA ---
for(i in 1:length(list.params)){
    key <- list.params[[i]]$key
    df <- filter_data(list.dat, r=list.params[[i]]$r, pval=list.params[[i]]$pval, correct_pval=list.params[[i]]$correct_pval)

    # WRITE OUTPUT ---
    file.output <- file.path(dir.tfpeaks2gene, paste("wcdt_tfpeaks2genes_correlated_", key, ".tsv", sep="") )
    write.table(df, file.output, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

    # COMPRESS ---
    cmd <- paste("gzip", file.output, sep=" ")
    system(cmd)
}





####################### GET distTSS #######################
### FUNCTION: filter_by_distTSS() ---
filter_by_distTSS <- function(annot, file.dat, threshold_kb){
    # LOAD DATA ---
    dat <- data.table::fread(file=file.dat, sep="\t", header=TRUE, nThread=60, data.table=FALSE, verbose=FALSE)

    # ADD distTSS ---
    df_annot <- subset(annot, annot$FeatureID %in% unique(dat$Peak))
    dat$distTSS <- NA
    for(i in 1:nrow(df_annot) ){
        index <- which( dat$Peak == df_annot$FeatureID[i] )
        dat$distTSS[index] <- dat_annot$distanceToTSS[i]
    }

    # FILTER BY distTSS ---
    dat$distTSS_kb <- round(dat$distTSS/1000, 2)
    dat <- dat[which(abs(dat$distTSS_kb) <= threshold_kb),]
    dat$distTSS_kb <- NULL

    return(dat)
}

### LOAD ATAC-SEQ READ COUNTS DATA ---
list.rds_atacseq_read_counts_mcrpc <- readRDS(file=file.rds_atacseq_read_counts_mcrpc)
annot <- list.rds_atacseq_read_counts_mcrpc$feature_annotation

### FILTER CORRELATION BY distTSS ---------------------------------------
file.dat <- file.path(dir.tfpeaks2gene, "wcdt_tfpeaks2genes_correlated_r_0p4_pval_0p05.tsv.gz")
df_500 <- filter_by_distTSS(annot, file.dat=, threshold_kb=500)
df_1000 <- filter_by_distTSS(annot, file.dat, threshold_kb=1000)

### WRITE OUTPUT ---
file.output_500 <- file.path(dir.tfpeaks2gene, "wcdt_tfpeaks2genes_correlated_r_0p4_pval_0p05_distTSS_500kb.tsv")
write.table(df_500, file.output_500, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
cmd_500 <- paste("gzip", file.output_500, sep=" ")
system(cmd_500)

### WRITE OUTPUT ---
file.output_1000 <- file.path(dir.tfpeaks2gene, "wcdt_tfpeaks2genes_correlated_r_0p4_pval_0p05_distTSS_1000kb.tsv")
write.table(df_1000, file.output_1000, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
cmd_1000 <- paste("gzip", file.output_1000, sep=" ")
system(cmd_1000)

