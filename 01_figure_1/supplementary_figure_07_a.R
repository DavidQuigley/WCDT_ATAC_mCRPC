###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### LOAD LIBRARIES ---
library("stringr")
library("data.table")
library("dplyr")
library("reshape2")
library("ggplot2")
library("gridExtra")

### DEFINE PATH ---
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

dir.bigwig <- file.path(dir.wrk, "secondary_data/wcdt_atacseq/bigwig") # REPLACE PATH
dir.data <- tempdir() # REPLACE PATH
dir.scripts <- tempdir() # REPLACE PATH
dir.logs <- tempdir() # REPLACE PATH

### DEFINE FILES ----
file.rds_atacseq_masterdata <- file.path(dir.reproduce_data, "wcdt_atacseq_masterdata.rds")
file.rds_rnaseq <- file.path(dir.reproduce_data, "wcdt_rnaseq_atacseq_gene_expression.rds")
file.utility_functions <- file.path(dir.reproduce_scripts, "00_prepare_data/00_01_utility_functions.R")

### CHROMOSOMES ---
chromosomes <- paste("chr", c(1:22,"X","Y"), sep="")

### LOAD FUNCTIONS ---
source(file.utility_functions)



### FUNCTION: get.SGEcmd()---
get.SGEcmd <- function(jobid, dir.logs){
    jobid <- jobid
    line01 <- "#!/bin/bash"
    line02 <- "#$ -pe smp 5"
    line03 <- "#$ -V"
    line04 <- "#$ -R y"
    line05 <- "#$ -l mem_free=10G"
    line06 <- paste("#$ -N", jobid, sep=" ")
    line07 <- paste("#$ -o", file.path(dir.logs, paste(jobid, ".out", sep="")), sep=" ")
    line08 <- paste("#$ -e", file.path(dir.logs, paste(jobid, ".err", sep="")), sep=" ")

    # COMPILE CMD ---
    cmd <- list(line01, line02, line03, line04, line05, line06, line07, line08)

    return(cmd)    
}



### FUNCTION: getComuteMatrix() ---
getComuteMatrix <- function(files.bigwig, file.bed, file.matrix, sampleids, n_cores=10){
    # 'computeMatrix' FUNCTION FROM 'deepTools'
    # https://deeptools.readthedocs.io/en/latest/index.html

    cmd <- paste("computeMatrix reference-point",
                "--scoreFileName", files.bigwig,
                "--regionsFileName", file.bed,
                "--outFileName", file.matrix,
                "--referencePoint", "TSS",
                "--beforeRegionStartLength", 3000,
                "--afterRegionStartLength", 3000,
                "--binSize", 100,
                "--samplesLabel", sampleids,
                "--missingDataAsZero",                
                "--verbose",
                "-p", n_cores,  
                sep=" ")
    
    return(cmd)
}



### FUNCTION: loadMatrix() ----
loadMatrix <- function(file.matrix){
    # LOAD MATRIX ---
    dat <- data.table::fread(file=file.matrix, sep="\t", header=FALSE, nThread=55, skip=1, data.table=FALSE, verbose=FALSE)
    rownames(dat) <- dat$V4
    dat$V1 <- dat$V2 <- dat$V3 <- dat$V4 <- dat$V5 <- dat$V6 <- NULL

    # ADD COLUMNS ---
    bins <- c(seq(-3000, -100, by=100), seq(100,3000, by=100))
    colnames(dat) <- bins

    #dat[is.na(dat)] <- 0

    return(dat)
}


### FUNCTION: getMeanProfilebyBED() ---
getMeanProfilebyBED <- function(dat, file.bed){
    # LOAD GENE BED DATA ---
    bed <- data.table::fread(file=file.bed, sep="\t", header=FALSE, nThread=1, data.table=FALSE, verbose=FALSE)

    # SUBSET MATRIX DATA BY GENES ---
    mat <- subset(dat, rownames(dat) %in% bed$V4)

    # GET MEAN PROFILE ---
    dat_mean <- apply(mat, 2, function(x) mean(x, na.rm=TRUE))

    return(dat_mean)
}

### FUNCTION: prepareData() ---
prepareData <- function(dat, dq){
    list.mat <- list()
    for(i in 1:nrow(dq)){
        q <- dq$Quantile[i]
        file.bed <- dq$File_bed[i]
        list.mat[[q]] <- getMeanProfilebyBED(dat=dat, file.bed=file.bed)
    }

    # AGGREGATE DATA ---
    mat <- do.call(rbind.data.frame, list.mat)
    rownames(mat) <- dq$Quantile
    colnames(mat) <- colnames(dat)
    
    # RESHAPE DATA ---
    dm <- reshape2::melt(as.matrix(mat))
    colnames(dm) <- c("Quantile","Position","Profile")

    # FACTORIZE ---
    dm$Quantile <- factor(dm$Quantile, levels=paste("Q", c(1:5), sep=""))

    return(dm)
}






######################################################################
### LOAD ATAC-SEQ MASTER DATA ---
list.rds_atacseq_masterdata <- readRDS(file=file.rds_atacseq_masterdata)
metadata <- list.rds_atacseq_masterdata$metadata_mcrpc

### DEFINE SAMPLE ---
sampleid <- "DTB-063-BL"




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


### ASSIGN QUANTILES ---
qexpr <- matrix(NA, nrow=nrow(expr), ncol=ncol(expr), dimnames=list(rownames(expr), colnames(expr)))

for(i in 1:ncol(expr)){
    qexpr[,i] <- dplyr::ntile(as.numeric(expr[,i]), 5)
}


### RE-ORDER ANNOTATION DATA ---
gannot <- gannot[match(rownames(expr), gannot$Gene),]


### EXTRACT GENE ANNOTATION DATA WITH QUANTILE INFO FOR THE SAMPLE ----
df <- gannot
df$ExprQuantile <- as.numeric(qexpr[,sampleid])
df <- df[order(df$ExprQuantile, decreasing=TRUE),]
df <- subset(df, select=c("Chrom","Start","End","Gene","ExprQuantile","Strand"))


### OUTPUT BED BY QUANTILE ---
list.files.bed <- list()
for(i in 1:5){
    quant <- i
    bed <- df[which(df$ExprQuantile == quant),]
    
    list.files.bed[[i]] <- file.path(dir.data, paste("genes_", sampleid, "_q_", quant, ".bed", sep=""))
    write.table(bed, list.files.bed[[i]], sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
}




##############################################################################################################################
### DEFINE INPUT FILES ---
file.bw <- file.path(dir.bigwig, paste(metadata$Sample_Name[which(metadata$Sample_ID == sampleid)], ".bw", sep=""))
files.bed <- paste(unlist(list.files.bed), collapse=" ")

### COMPUTE MATRIX ATAC-SEQ ---
file.matrix <- file.path(dir.data, "expr_quantiles_promoter_atacseq_matrix.tsv.gz")
cmd_1 <- getComuteMatrix(files.bigwig=file.bw, file.bed=files.bed, file.matrix, sampleids="DTB-063-BL", n_cores=40)

system(cmd_1)




##############################################################################################################################
### BED FILES ---
files.bed <- unlist(list.files.bed)
dq <- data.frame(Quantile=paste("Q", c(1:5), sep=""), File_bed=files.bed)
dq <- dq %>% dplyr::mutate_all(as.character)

### GET PROFILE ATAC-SEQ ---
dat_atacseq <- loadMatrix(file.matrix=file.matrix)
dm_atacseq <- prepareData(dat=dat_atacseq, dq=dq)


##############################################################################################################################
### PLOT ---
p_2 <- ggplot(data = dm_atacseq, aes(x=Position, y=Profile)) +       
                geom_line(aes(color=Quantile), size=0.6) +
                scale_x_continuous(breaks=c(-3000,-2000,-1000,0,1000,2000,3000), labels=c("-3","-2","-1","TSS","1","2","3")) +
                scale_color_manual(values=c("#fcbba1","#fc9272","#ef3b2c","#cb181d","#67000d")) +
                scale_y_continuous(breaks=seq(0,1.25, by=0.25)) + #ATAC-SEQ
                coord_cartesian(ylim=c(0,1.25)) + #ATAC-SEQ
                theme(
                    axis.text.x = element_text(size = 6, color="#000000"),
                    axis.text.y = element_text(size = 6, color="#000000"),
                    axis.title = element_text(size = 6, color="#000000"),
                    plot.title = element_text(size = 6, color="#000000", hjust=0.5),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    axis.ticks = element_line(size=0.2, color="#000000"),	
                    strip.text = element_text(size=6, color="#000000"),
                    strip.background = element_rect(fill="#FFFFFF", color="#FFFFFF"),
                    panel.background = element_rect(fill="#FFFFFF", color="#000000"),
                    legend.text = element_text(size = 6, color="#000000"),
                    legend.title = element_blank(),
                    legend.key.size = unit(0.3, "cm"),
                    legend.position = "none") +
                ylab("ATAC-seq Profile") +
                xlab("Distance From TSS") + 
                ggtitle("")



# WRITE PLOT ---
file.plot <- file.path(dir.reproduce_fig, "supplementary_figure_07_a.pdf")
pdf(file.plot, height=2, width=2)
    grid.arrange(p_2, nrow=1, ncol=1) 
dev.off()




