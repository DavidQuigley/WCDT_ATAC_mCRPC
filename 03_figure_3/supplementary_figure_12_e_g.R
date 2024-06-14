###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### LOAD LIBRARIES ---
library("stringr")
library("data.table")
library("reshape2")
library("ggplot2")
library("gridExtra")

### DEFINE PATH ---
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/analysis/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

### DEFINE FILES ----
file.rds_atacseq_masterdata <- file.path(dir.reproduce_data, "wcdt_atacseq_masterdata.rds")
file.rds_atacseq_read_counts_mcrpc <- file.path(dir.reproduce_data, "atacseq_read_counts_mcrpc.rds")
file.rds_rnaseq <- file.path(dir.reproduce_data, "wcdt_rnaseq_atacseq_gene_expression.rds")
file.kwtest <- file.path(dir.reproduce_data, "wcdt_atacseq_mcrpc_subtypes_KruskalWallisTest_summary.tsv.gz")
file.utility_functions <- file.path(dir.reproduce_scripts, "00_prepare_data/00_01_utility_functions.R")

file.gene_annotation <- file.path(dir.wrk, "reference/annotation_data/gene_annotation/gencode_v28/gencode.v28.annotation.gtf.gz")

######################################################################
### SUBTYPES ---
subtypes <- c("AR+_NE-","ARL_NE-","AR+_NE+","AR-_NE+","AR-_NE-")

### COLOR ---
cp_subtypes <- c("#e31a1c","#fb9a99","#ff7f00","#33a02c","#1f78b4")

### LOAD GENE FEATURE DATA --
genes <- c("AR","NKX3-1","GPR37L1")

### LOAD FUNCTIONS ---
source(file.utility_functions)



############################################################################################################################################
### FUNCTION: getPlotRange() ---
getPlotRange <- function(annot.gene, gene, upstream, downstream){
    # GET GENE POSITION ---
    plot_chr <- annot.gene$Chrom[which(annot.gene$Gene == gene)] 
    gene_start <- annot.gene$Start[which(annot.gene$Gene == gene)] 
    gene_end <- annot.gene$End[which(annot.gene$Gene == gene)]
    plot_region_start <- gene_start - upstream   #50000 bp
    plot_region_end <- gene_end + downstream     #10000 bp

    # RANGE ---
    gr_range_limit <- GenomicRanges::GRanges(
                                seqnames = Rle(plot_chr),
                                ranges = IRanges(start=plot_region_start, end=plot_region_end),
                                strand = Rle(strand("+")))    
    
    return(gr_range_limit)
}

### FUNCTION: parseDiffPeaks() ---
parseDiffPeaks <- function(gr_diffpeaks){
    df_diffpeaks <- as.data.frame(gr_diffpeaks)
    df_diffpeaks$FeatureID <- apply(df_diffpeaks, 1, function(x) paste(x[1], as.numeric(x[2]), as.numeric(x[3]), sep="_"))
    rownames(df_diffpeaks) <- NULL
    return(df_diffpeaks)
}

### FUNCTION: getCorr() ---
getCorr <- function(x, y){
    correlation.test <- cor.test(x, y, method="spearman", na.action = "na.exclude")
	r <- round(as.numeric(correlation.test$estimate),2)
    return(r)
}	



### FUNCTION: prepareData() ---
prepareData <- function(mat_atacseq, expr, feature_id, subtype, gene){
    # PREPARE DATA ---
    dm <- data.frame(SampleID=colnames(mat_atacseq), 
                    FeatureID=feature_id, Gene=gene, Subtype=subtype,
                    ATACseq=as.numeric(mat_atacseq[feature_id,]), 
                    RNAseq=as.numeric(expr[gene,]))

    return(dm)
}


### FUNCTION: getPlot() ---
getPlot <- function(dm, feature_id, subtypes, colors){
    # FACTORIZE ---
    dm$Subtype <- factor(dm$Subtype, levels=subtypes)

    # PLOT ---
    p <- ggplot(dm, aes(x=ATACseq, y=RNAseq)) + 
            stat_smooth(method="lm", geom = "smooth", formula = y ~ x, position = "identity", size=0.5, fullrange = FALSE, se=TRUE, na.rm=TRUE) +	
            geom_point(aes(fill=Subtype, color=Subtype), shape = 21, stroke=0.5, size=0.8, alpha=0.7, na.rm=TRUE) +
            #geom_point(fill="#e41a1c", color="#e41a1c", shape = 21, stroke=0.5, size=0.8, alpha=0.7, na.rm=TRUE) +
            #coord_cartesian(xlim=c(20000, 160000), ylim=c(0, 0.3)) +
            #scale_x_continuous(breaks=seq(20000,160000, by=20000), labels=seq(20000,160000, by=20000)) +
            #scale_y_continuous(breaks=seq(0,0.3, by=0.05), labels=seq(0,0.3, by=0.05) ) +
            #annotate("text", x = 120000, y = 0.025, label = paste("R = ", r, sep=""), size=2.5) +	
            scale_fill_manual(values=colors) +
            scale_color_manual(values=colors) +
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
            ylab("RNA-seq") +
            xlab("ATAC-seq") + 
            ggtitle(feature_id) 

    return(p)
}


### FUNCTION: getCorrPlot() ---
getCorrPlot <- function(metadata, dat_atacseq, expr, gr_diffpeaks, annot.gene, subtypes, cp_subtypes, gene, upstream, downstream){
    # GET PLOT RANGE ---
    gr_range_limit <- getPlotRange(annot.gene, gene=gene, upstream, downstream)

    # DIFF ATAC-SEQ PEAKS ---
    gr_diffpeaks_plot <- IRanges::subsetByOverlaps(x=gr_diffpeaks, ranges=gr_range_limit, type="any")
    df_diffpeaks <- parseDiffPeaks(gr_diffpeaks_plot)

    # EXTRACT FEATURES FROM ATAC-SEQ COUNTS MATRIX ---
    mat_atacseq <- subset(dat_atacseq, rownames(dat_atacseq) %in% df_diffpeaks$FeatureID)

    # GET CORRELATION ----
    list.df <- list()
    for(i in 1:nrow(df_diffpeaks)){
        feature_id <- df_diffpeaks$FeatureID[i]

        x <- as.numeric(mat_atacseq[feature_id,])
        y <- as.numeric(expr[gene,])
        r <- getCorr(x, y)

        list.df[[feature_id]] <- data.frame(feature_id=feature_id, gene=gene, r=r)
    }

    # AGGREGATE DATA ---
    df <- do.call(rbind.data.frame, list.df)
    rownames(df) <- NULL


    # PREPARE DATA FOR EACH FEATURE ---
    list.dm <- list()
    for(i in 1:nrow(df_diffpeaks)){
        feature_id <- df_diffpeaks$FeatureID[i]
        list.dm[[feature_id]] <- prepareData(mat_atacseq, expr, feature_id=feature_id, subtype=metadata$mCRPC_Subtype, gene=gene)
    }

    # GET PLOT ---
    list.plot <- list()
    for(i in 1:nrow(df_diffpeaks)){
        feature_id <- df_diffpeaks$FeatureID[i]
        list.plot[[feature_id]] <- getPlot(dm=list.dm[[feature_id]], feature_id=feature_id, subtypes=subtypes, colors=cp_subtypes)
    }

    # COMPILE OUTPUT ---
    list.output <- list(datcorr=df, dat=list.dm, plot=list.plot)

    return(list.output)
}

############################################################################################################################################



######################################################################
### LOAD ATAC-SEQ MASTER DATA ---
list.rds_atacseq_masterdata <- readRDS(file=file.rds_atacseq_masterdata)
metadata <- list.rds_atacseq_masterdata$metadata_mcrpc

### LOAD ATAC-SEQ READ COUNTS DATA ---
list.rds_atacseq_read_counts_mcrpc <- readRDS(file=file.rds_atacseq_read_counts_mcrpc)
dat_atacseq <- list.rds_atacseq_read_counts_mcrpc$feature_counts_norm
rownames(dat_atacseq) <- dat_atacseq$feature_id
dat_atacseq$feature_id <- NULL
colnames(dat_atacseq) <- metadata$Sample_ID

### LOAD ATAC-SEQ PEAK ANNOTATION DATA ---
annot <- list.rds_atacseq_read_counts_mcrpc$feature_annotation

### GET REF ANNOTATION DATA FOR GENE AND EXON ---
annot.gene <- parseGTF(file.gtf=file.gene_annotation, feature.type="gene")



### LOAD KRUSKAL-WALLIS TEST DATA ---
dat_kwtest <- data.table::fread(file=file.kwtest, sep="\t", header=TRUE, nThread=50, data.table=FALSE, verbose=FALSE)
colnames(dat_kwtest)[1] <- "FeatureID"
dat_kwtest <- merge(dat_kwtest, annot, by="FeatureID")
dat_kwtest <- dat_kwtest[order(dat_kwtest$pvalue, decreasing=FALSE),]
dat_kwtest <- subset(dat_kwtest, dat_kwtest$pvalue <= 0.001)

dat_kwtest$chr <- unlist(lapply(stringr::str_split(dat_kwtest$FeatureID, "_"), function(x) x[1]))
dat_kwtest$start <- as.numeric(unlist(lapply(stringr::str_split(dat_kwtest$FeatureID, "_"), function(x) x[2])))
dat_kwtest$end <- as.numeric(unlist(lapply(stringr::str_split(dat_kwtest$FeatureID, "_"), function(x) x[3])))

diffpeaks_bed <- subset(dat_kwtest, select=c("chr","start","end"))
gr_diffpeaks <- GenomicRanges::makeGRangesFromDataFrame(df=diffpeaks_bed, keep.extra.columns=TRUE)





################## RNA-SEQ DATA ################################################################################################

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




### GET PLOTS ----
list.genes1 <- getCorrPlot(metadata, dat_atacseq, expr, gr_diffpeaks, annot.gene, subtypes, cp_subtypes, gene="AR", upstream=800000, downstream=-100000)
list.genes2 <- getCorrPlot(metadata, dat_atacseq, expr, gr_diffpeaks, annot.gene, subtypes, cp_subtypes, gene="NKX3-1", upstream=12000, downstream=5000)
list.genes3 <- getCorrPlot(metadata, dat_atacseq, expr, gr_diffpeaks, annot.gene, subtypes, cp_subtypes, gene="GPR37L1", upstream=60000, downstream=10000)


#> list.genes1$datcorr
#              feature_id gene    r
#1 chrX_67019495_67019638   AR 0.75
#2 chrX_66899614_66900940   AR 0.86
#3 chrX_67622403_67623063   AR 0.83
#4 chrX_67039188_67039632   AR 0.83
#5 chrX_67008634_67009006   AR 0.76
#6 chrX_67612181_67612655   AR 0.80

#> list.genes2$datcorr
#              feature_id   gene    r
#1 chr8_23677512_23679924 NKX3-1 0.74
#2 chr8_23680335_23680906 NKX3-1 0.54
#3 chr8_23674041_23674454 NKX3-1 0.60
#4 chr8_23671732_23672422 NKX3-1 0.67

#> list.genes3$datcorr
#                feature_id    gene    r
#1 chr1_202108098_202108759 GPR37L1 0.24
#2 chr1_202105754_202106079 GPR37L1 0.11
#3 chr1_202107120_202107987 GPR37L1 0.23
#4 chr1_202112630_202113793 GPR37L1 0.27
#5 chr1_202106317_202107052 GPR37L1 0.39
#6 chr1_202101716_202101941 GPR37L1 0.38
#7 chr1_202084763_202086340 GPR37L1 0.36


### GENERATE PLOT ---
file.plot <- file.path(dir.reproduce_fig, "supplementary_figure_12_e_g.pdf")
pdf(file.plot, width=4, height=1.5)
    grid.arrange(list.genes1$plot[["chrX_66899614_66900940"]], nrow=1, ncol=3)
    grid.arrange(list.genes2$plot[["chr8_23674041_23674454"]], list.genes2$plot[["chr8_23671732_23672422"]], nrow=1, ncol=3)
    grid.arrange(list.genes3$plot[["chr1_202106317_202107052"]], list.genes3$plot[["chr1_202101716_202101941"]], list.genes3$plot[["chr1_202084763_202086340"]], nrow=1, ncol=3)
dev.off()
