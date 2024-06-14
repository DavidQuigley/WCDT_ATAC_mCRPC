###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### LOAD LIBRARIES ---
library("stringr")
library("data.table")
library("FNN")
library("ggplot2")
library("gridExtra")

### DEFINE PATH ---
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/analysis/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

### DEFINE FILES ----
file.rds_atacseq_masterdata <- file.path(dir.reproduce_data, "wcdt_atacseq_masterdata.rds")
file.rds_atacseq_read_counts_combined <- file.path(dir.reproduce_data, "atacseq_read_counts_combined.rds")
file.rds_rnaseq <- file.path(dir.reproduce_data, "tcga_wcdt_rnaseq_gene_expression.rds")
file.utility_functions <- file.path(dir.reproduce_scripts, "00_prepare_data/00_01_utility_functions.R")

### LOAD FUNCTIONS ---
source(file.utility_functions)


### FUNCTION: get.wilcox.rank.test() ----
get.wilcox.rank.test <- function(dat, class1, class2){
	cat("COMPUTE WILCOX RANK TEST ...", "\n", sep="")
	#dat <- dat[-which(is.na(dat)),]
	
	#Get class indices
	index1 <- which(colnames(dat) %in% class1)
	index2 <- which(colnames(dat) %in% class2)
	
	# Check if the data points are constant ----
	exp <- dat[,c(index1, index2)]
	len <- apply(exp, 1, function(x) length(unique(as.numeric(x))))
	del.index <- which(len == 1)

	if(length(del.index) != 0){
		dat <- dat[-del.index,]
	}

	# Check if each group has at least 2 samples represented ---
	length_not_na.A <- as.numeric(apply(dat[,index1], 1, function(x) length(which(!is.na(as.numeric(x))))))
	length_not_na.B <- as.numeric(apply(dat[,index2], 1, function(x) length(which(!is.na(as.numeric(x))))))
	del.index_A <- which(length_not_na.A < 2)
	del.index_B <- which(length_not_na.B < 2)

	del.index_AB <- unique(c(del.index_A, del.index_B))
	
	if(length(del.index_AB) != 0){
		dat <- dat[-del.index_AB,]
	}
	
	# Compute Median and StdDev for the Sample Class
	median1 <- apply(as.matrix(dat)[,index1], 1, median, na.rm=TRUE)
	median2 <- apply(as.matrix(dat)[,index2], 1, median, na.rm=TRUE)
	cat("MEDIAN COMPUTED ...", "\n", sep="")
	
	stdev1 <- apply(as.matrix(dat)[,index1], 1, sd, na.rm=TRUE)
	stdev2 <- apply(as.matrix(dat)[,index2], 1, sd, na.rm=TRUE)
	cat("STANDARD DEVIATION COMPUTED ...", "\n", sep="")
	
	# Compute FoldChange
	foldchange <- ifelse(median2 > median1, 2^median2/2^median1, -2^median1/2^median2)
	cat("FOLD CHANGE COMPUTED ...", "\n", sep="")
	
	cat("PERFORMING WILCOX RANK TEST ...", "\n", sep="")
	# Wilcox.rank.test
	func.wilcoxtest <- function(x, index1, index2){
		wilcoxtest <- suppressWarnings(wilcox.test(x[index1], x[index2], alternative="two.sided", na.action=na.omit))
		return(wilcoxtest)
	}

	wilcoxtest.list <- suppressWarnings(apply(as.matrix(dat), 1, function(x) func.wilcoxtest(x, index1, index2)))
	cat("WILCOX RANK TEST COMPUTED ...", "\n", sep="")
	
	pvalue <- suppressWarnings(do.call(c,lapply(wilcoxtest.list,function(x) x$p.value)))

	index0 <- which(pvalue == 0)
	if(length(index0) != 0){
		pvalue[index0] <- .Machine$double.eps #smallest value
	}

	fdr <- suppressWarnings(p.adjust(pvalue, method = "BH", n = length(pvalue)))
	cat("P-VALUE COMPUTED ...", "\n", sep="")

	dat.summary <- data.frame(Gene=rownames(dat), 
								MedianA=median1, StdevA=stdev1, 
								MedianB=median2, StdevB=stdev2, 
								FoldChange=foldchange, 
								pvalue=pvalue, fdr=fdr)
	dat.summary$Gene <- as.character(dat.summary$Gene)
	dat.summary <- dat.summary[order(dat.summary$pvalue, decreasing=FALSE),]
    rownames(dat.summary) <- NULL

	cat("DATA SUMMARY COMPUTED ...", "\n", sep="")
		
	return(dat.summary)
}

### FUNCTION: prepareData() ---
prepareData <- function(diff_atac, diff_rna){
    # GET COMMON GENES ----
    genes <- intersect(unique(diff_atac$Gene), unique(diff_rna$Gene))

    # EXTRACT FODCHANGE FOR EACH GENE ---
    list.df <- list()
    for(i in 1:length(genes)){
        gene <- genes[i]

        # GET FOLDCHANGE VALUES ---
        val_rna <- diff_rna$log2FoldChange[which(diff_rna$Gene == gene)]
        vals_atac <- diff_atac$log2FoldChange[which(diff_atac$Gene == gene)]

        # GET NEAREST ATAC-SEQ FC VALUE TO THAT OF RNA-SEQ ---
        if(length(vals_atac) > 1){    
            index_knn <- FNN::get.knnx(data=vals_atac, query=val_rna, k=1, algorithm="kd_tree")$nn.index
            val_atac <- vals_atac[index_knn]
        }else{
            val_atac <- vals_atac
        }

        # PREPARE DATA ---
        list.df[[gene]] <- data.frame(Gene=gene, log2fc_atac=val_atac, log2fc_rna=val_rna)
    }

    # AGGREGATE DATA ---
    df <- do.call(rbind.data.frame, list.df)
    rownames(df) <- NULL

    # ADD SIGN ---
    df$sign_atac <- sign(df$log2fc_atac)
    df$sign_rna <- sign(df$log2fc_rna)

    # ADD CORR ---
    df$corr <- ifelse(( (df$sign_atac == 1) & (df$sign_rna == 1) ) | ( (df$sign_atac == -1) & (df$sign_rna == -1) ), 1, 0)

    return(df)
}

### FUNCTION: getCorr() ---
getCorr <- function(x, y){
    correlation.test <- cor.test(x, y, method="spearman", na.action = "na.exclude")
	r <- round(as.numeric(correlation.test$estimate),2)
    return(r)
}

### FUNCTION: getPlot() ---
getPlot <- function(df, feature_type){
    p <- ggplot(df, aes(x=log2fc_atac, y=log2fc_rna)) +
    		geom_hline(yintercept = 0, color="#000000", linetype=4, alpha=0.7) +
			geom_vline(xintercept = 0, color="#000000", linetype=4, alpha=0.7) +		 
            geom_hex(fill="#e41a1c", color="#e41a1c", stat="binhex", position="identity", bins=300, alpha=0.3) +
            stat_smooth(method="lm", geom = "smooth", formula = y ~ x, position = "identity", size=0.5, fullrange = FALSE, se=TRUE, na.rm=TRUE) +	
            #geom_point(fill="#e41a1c", color="#e41a1c", shape = 21, stroke=0.5, size=0.5, alpha=0.5, na.rm=TRUE) +
            #coord_cartesian(xlim=c(-8, 8), ylim=c(-8, 8)) +
            #scale_x_continuous(breaks=seq(-8,8, by=2), labels=seq(-8,8, by=2)) +
            #scale_y_continuous(breaks=seq(-8,8, by=2), labels=seq(-8,8, by=2) ) +
            theme(
                axis.text.x = element_text(size = 5, color="#000000"),
                axis.text.y = element_text(size = 5, color="#000000"),
                axis.title = element_text(size = 5, color="#000000"),
                plot.title = element_text(size = 5, color="#000000", hjust=0.5),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.ticks = element_line(size=0.2, color="#000000"), 
                strip.text = element_text(size=5, color="#000000"),
                strip.background = element_rect(fill="#FFFFFF", color="#FFFFFF"),
                panel.background = element_rect(fill="#FFFFFF", color="#000000"),
                legend.text = element_text(size = 5, color="#000000"),
                legend.title = element_blank(),
                legend.key.size = unit(0.3, "cm"),
                legend.position = "none") +
            ylab("log2FC (RNA-seq)") +
            xlab("log2FC (ATAC-seq)") + 
            ggtitle(feature_type) 

    return(p)                
}


##############################################################################################################################
### LOAD ATAC-SEQ MASTER DATA ---
list.rds_atacseq_masterdata <- readRDS(file=file.rds_atacseq_masterdata)

### METADATA: WCDT mCRPC ----
metadata <- list.rds_atacseq_masterdata$metadata_mcrpc
metadata <- subset(metadata, metadata$mCRPC_Subtype != "AR-_NE+")

### METADATA: COMBINED ---
metadata_combo <- list.rds_atacseq_masterdata$metadata_combined
metadata_combo <- subset(metadata_combo, metadata_combo$Dataset %in% c("TCGA","WCDT"))
metadata_combo <- subset(metadata_combo, metadata_combo$Phenotype %in% c("PCa","mCRPC-Adeno"))
metadata_combo <- subset(metadata_combo, select=c("SampleName","Phenotype","Dataset"))
metadata_combo <- rbind( data.frame(SampleName="TCGA-EJ-A46G-01A", Phenotype="PCa", Dataset="TCGA"), metadata_combo )

#> all(metadata$Sample_Name ==  metadata_combo$SampleName[which(metadata_combo$Dataset == "WCDT")])
#[1] TRUE

metadata_combo$SampleID <- NA
metadata_combo$SampleID[which(metadata_combo$Dataset == "TCGA")] <- metadata_combo$SampleName[which(metadata_combo$Dataset == "TCGA")]
metadata_combo$SampleID[which(metadata_combo$Dataset == "WCDT")] <- metadata$Sample_ID
metadata_combo <- subset(metadata_combo, select=c("SampleID","SampleName","Phenotype","Dataset"))
metadata_combo$Phenotype <- stringr::str_replace_all(metadata_combo$Phenotype, "mCRPC-Adeno", "mCRPC")


### GET SAMPLEIDS ---
sampleids_pca <- metadata_combo$SampleID[which(metadata_combo$Phenotype == "PCa")]
sampleids_mcrpc <- metadata_combo$SampleID[which(metadata_combo$Phenotype == "mCRPC")]

#"TCGA-EJ-A46G-01A"


##############################################################################################################################


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
    expr <- subset(expr, rownames(expr) %in% gannot$EnsemblID)
    rownames(expr) <- gannot$Gene
cat(format(Sys.time(), "%b %d %X"), "DONE! ", "\n", sep=" ")

### RE-ARRANGE EXPRESSION DATA ---
expr <- expr[order(rownames(expr), decreasing=FALSE),]
expr <- subset(expr, select=metadata_combo$SampleID)

### LOG2 TRANSFORM [log2(TPM + 1)]---
#expr <- log2(expr + 1)

########## -------------------------------------------------------------------------------------------------------------------


### GET DIFFERENTIAL RNA-SEQ GENE-EXPRESSION ---
df_rnaseq <- get.wilcox.rank.test(dat=expr, class1=sampleids_pca, class2=sampleids_mcrpc)
    

### EXTRACT DIFFEXPR GENES ---
df_rnaseq_filtered <- subset(df_rnaseq, df_rnaseq$fdr <= 0.01)
df_rnaseq_filtered <- subset(df_rnaseq_filtered, (df_rnaseq_filtered$FoldChange <= -2) | (df_rnaseq_filtered$FoldChange >= 2) )
df_rnaseq_filtered$sign <- sign(df_rnaseq_filtered$FoldChange)
df_rnaseq_filtered$log2FoldChange <- apply(df_rnaseq_filtered, 1, function(x) log2(abs(as.numeric(x[6]))) * as.numeric(x[9]))  

#> dim(df_rnaseq_filtered)
#[1] 6388   10




######## ATAC-SEQ DATA ############################################################################################################
### LOAD ATAC-SEQ READ COUNTS DATA ---
list.rds_atacseq_read_counts_combined <- readRDS(file=file.rds_atacseq_read_counts_combined)
dat_atacseq <- list.rds_atacseq_read_counts_combined$feature_counts_raw
rownames(dat_atacseq) <- dat_atacseq$feature_id
dat_atacseq$feature_id <- NULL

### METADATA ---
metadata_combo <- metadata_combo[-which(metadata_combo$SampleID == "TCGA-EJ-A46G-01A"),]

### TRIM DATA BY SAMPLEID ---
dat_atacseq <- subset(dat_atacseq, select=metadata_combo$SampleName )
colnames(dat_atacseq) <- metadata_combo$SampleID



### GET DESEQ OBJECT ---
dds <- DESeq2::DESeqDataSetFromMatrix(countData=dat_atacseq, colData=metadata_combo, design= ~ Phenotype)
dds <- DESeq2::DESeq(dds)

### FUNCTION: getDiffExpression() ---
getDiffExpression <- function(obj, groupA, groupB){
    item <- c("Phenotype", groupA, groupB)
    res <- DESeq2::results(obj, contrast=item, independentFiltering=FALSE, pAdjustMethod="BH")

    index0 <- which(res$pvalue == 0)
    if(length(index0) != 0){
        res$pvalue[index0] <- .Machine$double.eps #smallest value
        res$padj <- p.adjust(res$pvalue, method = "BH", n = length(res$pvalue))
    }

    return(res)
}


### DIFF. PEAKS: PCa vs Benign ---
diffPeak <- getDiffExpression(obj=dds, groupA="mCRPC", groupB="PCa")    
diffPeak <- as.data.frame(diffPeak)
diffPeak$FeatureID <- rownames(diffPeak)
rownames(diffPeak) <- NULL

### MERGE DIFFPEAKS AND ANNOT DATA ---
df_atacseq <- merge(diffPeak, list.rds_atacseq_read_counts_combined$feature_annotation, by="FeatureID")

### SORT DATA BY FDR VALUE ----
df_atacseq <- df_atacseq[order(df_atacseq$padj, decreasing=FALSE),]
rownames(df_atacseq) <- NULL



### FUNCTION: filterData() ---
filterData <- function(dat, pvalue.threshold=0.01, logfc.threshold=1){
    # DIFFERENTIAL ACCESSIBILITY ---
    d <- subset(dat, dat$padj <= pvalue.threshold)
    d <- subset(d, (d$log2FoldChange >= logfc.threshold) | (d$log2FoldChange <= -logfc.threshold) )

    return(d)
}

### FILTER DIFF-PEAK DATA ---
df_atacseq_filtered <- filterData(dat=df_atacseq, pvalue.threshold=0.01, logfc.threshold=1)
df_atacseq_prot <- subset(df_atacseq_filtered, df_atacseq_filtered$GeneType == "protein_coding")

#> dim(df_atacseq_prot)
#[1] 16599    16

#> length(unique(df_atacseq_prot$Gene))
#[1] 7343




######################################################################################################################################
### GET UNION OF DIFF PEAKS AND EXPR GENES ---
genes_diff_atac_u_rna <- unique(c(unique(df_atacseq_prot$Gene), unique(df_rnaseq_filtered$Gene)))

#> length(genes_diff_atac_u_rna)
#[1] 11221


### SUBSET ATAC AND RNA DATA BY DIFF GENES ---
d_atac <- subset(df_atacseq, df_atacseq$Gene %in% genes_diff_atac_u_rna)

d_rna <- subset(df_rnaseq, df_rnaseq$Gene %in% genes_diff_atac_u_rna)
d_rna$sign <- sign(d_rna$FoldChange)
d_rna$log2FoldChange <- apply(d_rna, 1, function(x) log2(abs(as.numeric(x[6]))) * as.numeric(x[9]))  

### PREPARE DATA ---
df_union_genes <- prepareData(diff_atac=d_atac, diff_rna=d_rna)

### GET SPEARMAN'S CORRELATION ---
r <- getCorr(x=df_union_genes$log2fc_atac, y=df_union_genes$log2fc_rna)

### GET PLOT ---
p_union_genes <- getPlot(df=df_union_genes, feature_type="ATAC-seq RNA-seq Diff. union genes")

### GENERATE PLOT ---
file.plot <- file.path(dir.reproduce_fig, "supplementary_figure_07_c.pdf")
pdf(file.plot, width=2, height=2)
    grid.arrange(p_union_genes, nrow=1, ncol=1)
dev.off()


