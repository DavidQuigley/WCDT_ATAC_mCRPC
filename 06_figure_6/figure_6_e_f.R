###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### LOAD LIBRARIES ---
library("stringr")
library("data.table")
library("GenomicRanges")
library("reshape2")
library("ggplot2")
library("gridExtra")

### DEFINE PATH ---
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/analysis/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

dir.footprints <- file.path(dir.reproduce_data, "tf_footprints/mcrpc_samplewise")

### DEFINE FILES ---
file.rds_atacseq_masterdata <- file.path(dir.reproduce_data, "wcdt_atacseq_masterdata.rds")
file.rds_rnaseq <- file.path(dir.reproduce_data, "wcdt_rnaseq_atacseq_gene_expression.rds")
file.rds_genelist <- file.path(dir.reproduce_data, "genelist.rds")
file.rds_footprints_bound <- file.path(dir.footprints, "data_footprints/mcrpc_samplewise_motif_position_bound.rds")
file.rds_footprints_unbound <- file.path(dir.footprints, "data_footprints/mcrpc_samplewise_motif_position_unbound.rds")
file.annotation_gtf <- file.path(dir.wrk, "reference/annotation_data/gene_annotation/gencode_v28/gencode.v28.annotation.gtf.gz")
file.utility_functions <- file.path(dir.reproduce_scripts, "00_prepare_data/00_01_utility_functions.R")

### LOAD FUNCTIONS ---
source(file.utility_functions)

### SUBTYPES ---
subtypes <- c("ARpNEn","ARlNEn","ARpNEp","ARnNEp","ARnNEn")

### MOTIF ---
motif_prefix_znf <- "MA0528.1_ZNF263"
motif_prefix_myc <- "MA0147.3_MYC"



################################################################################################################################################
### FUNCTION: getFootprintRanges() ---
getFootprintRanges <- function(feature_ids){
    # PREPARE DATA ---
    df <- data.frame(FeatureID=feature_ids)
    
    # GET GENOMIC COORDINATES ---
    df$chr <- unlist(lapply(stringr::str_split(df$FeatureID, "_"), function(x) x[1]))
    df$start <- as.numeric(unlist(lapply(stringr::str_split(df$FeatureID, "_"), function(x) x[2])))
    df$end <- as.numeric(unlist(lapply(stringr::str_split(df$FeatureID, "_"), function(x) x[3])))

    # GET GENOMIC RANGES OBJECT ---
    gr <- GenomicRanges::makeGRangesFromDataFrame(df, keep.extra.columns=TRUE)

    return(gr)
}


### FUNCTION: annotPeaksbyRange() ---
annotPeaksbyRange <- function(gr_peaks, genecode.txdb){
    # ANNOTATE PEAKS ---
    annot_peaks <- ChIPseeker::annotatePeak(peak=gr_peaks, TxDb=genecode.txdb, verbose=TRUE)
    annot_peaks <- as.data.frame(annot_peaks)

    # ADD FEATURES ---
    annot_peaks$Feature <- ""
    annot_peaks$Feature[which(annot_peaks$annotation == "Promoter (<=1kb)")] <- "Promoter"
    annot_peaks$Feature[which(annot_peaks$annotation == "Promoter (1-2kb)")] <- "Promoter"
    annot_peaks$Feature[which(annot_peaks$annotation == "Promoter (2-3kb)")] <- "Promoter"
    annot_peaks$Feature[which(annot_peaks$annotation == "5' UTR")] <- "5UTR"
    annot_peaks$Feature[which(annot_peaks$annotation == "3' UTR")] <- "3UTR"

    exon <- annot_peaks$annotation[which(stringr::str_detect(annot_peaks$annotation, "Exon") == TRUE)]
    annot_peaks$Feature[which(annot_peaks$annotation %in% exon[which(str_detect(exon, "exon 1 of") == TRUE)])] <- "Exon"
    annot_peaks$Feature[which(annot_peaks$annotation %in% exon[which(str_detect(exon, "exon 1 of") == FALSE)])] <- "Exon"

    intn <- annot_peaks$annotation[which(stringr::str_detect(annot_peaks$annotation, "Intron") == TRUE)]
    annot_peaks$Feature[which(annot_peaks$annotation %in% intn[which(str_detect(intn, "intron 1 of") == TRUE)])] <- "Intron"
    annot_peaks$Feature[which(annot_peaks$annotation %in% intn[which(str_detect(intn, "intron 1 of") == FALSE)])] <- "Intron"

    annot_peaks$Feature[which(stringr::str_detect(annot_peaks$annotation, "Downstream") == TRUE)] <- "Downstream"
    annot_peaks$Feature[which(annot_peaks$annotation == "Distal Intergenic")] <- "Intergenic"

    # ADD FEATURE ID ---
    annot_peaks$FeatureID <- apply(annot_peaks, 1, function(x) paste(x[1], as.numeric(x[2]), as.numeric(x[3]), sep="_"))

    # ARRANGE COLUMNS ---
    items <- c("FeatureID","geneId","transcriptId","Feature","distanceToTSS")
    annot_peaks <- subset(annot_peaks, select=items)

    return(annot_peaks)
}

### FUNCTION: addGenes() ---
addGenes <- function(annot, dat){
    pb <- txtProgressBar(min = 0, max = nrow(dat), style = 3)
    dat$Gene <- ""
    dat$GeneType <- ""
    n <- nrow(dat)

    # ADD GENE SYBBOLS ---
    for(i in 1:n){
        index <- which(annot$EnsemblTranscriptID == dat$transcriptId[i])
        if(length(index) != 0){
            dat$Gene[i] <- annot$Gene[index]
            dat$GeneType[i] <- annot$GeneType[index]

            # UPDATE PROGRESS BAR ---
            setTxtProgressBar(pb, i)
        }
    }

    return(dat)
}



### FUNCTION: getGeneLevelMatrix() ---   
getGeneLevelMatrix <- function(list.footprints_bound_motif, dat){
    # CREATE GENE LEVEL MATRIX FOR BOUND/UNBOUND IN PROMOTERS ---
    genes <- sort(unique(dat$Gene))
    mat_genes <- matrix(0, nrow=length(genes), ncol=length(list.footprints_bound_motif), dimnames=list(genes, names(list.footprints_bound_motif)))

    ### LABEL BOUND GENES IN THE MATRIX ----
    for(i in 1:length(list.footprints_bound_motif)){
        sampleid <- names(list.footprints_bound_motif)[i]
        df <- list.footprints_bound_motif[[sampleid]]
        f_ids <- unique(df$motif_position)

        genes_bound <- unique(dat$Gene[which(dat$FeatureID %in% f_ids)])

        index_bound <- which(rownames(mat_genes) %in% genes_bound)
        mat_genes[index_bound,i] <- 1

        cat("PROCESSED:", sampleid, "\n", sep=" ")
    }

    return(mat_genes)
}



### FUNCTION: prepareMatrixData() ---
prepareMatrixData <- function(list.footprints_bound, list.footprints_unbound, metadata, annot, genecode.txdb, motif_prefix){
    # EXTRACT MOTIF FOOTPRINTS ----
    list.footprints_bound_motif <- list.footprints_bound[[motif_prefix]]
    list.footprints_unbound_motif <- list.footprints_unbound[[motif_prefix]]

    # ADD SAMPLEIDS ---
    names(list.footprints_bound_motif) <- metadata$Sample_ID
    names(list.footprints_unbound_motif) <- metadata$Sample_ID

    # GET FEATUREIDS ---
    feature_ids_bound <- unique(unlist(lapply(list.footprints_bound_motif, function(x) x$motif_position)))
    feature_ids_unbound <- unique(unlist(lapply(list.footprints_unbound_motif, function(x) x$motif_position)))

    # GET UNIQUE FEATUREIDS ---
    feature_ids <- unique(c(feature_ids_bound, feature_ids_unbound))
    gr_footprints <- getFootprintRanges(feature_ids)

    # GET FOOTPRINT SITE ANNOTATION -----
    dat_footprint_annot <- annotPeaksbyRange(gr_peaks=gr_footprints, genecode.txdb)

    # ADD GENE INFO ---
    dat_footprint_annot <- addGenes(annot, dat=dat_footprint_annot)
    dat_footprint_annot <- subset(dat_footprint_annot, dat_footprint_annot$GeneType == "protein_coding")

    # GET GENE LEVEL MATRIX ---
    mat_gene <- getGeneLevelMatrix(list.footprints_bound_motif, dat=dat_footprint_annot)

    return(mat_gene)
}

### FUNCTION: prepareDataStats() ---
prepareDataStats <- function(df, dat_expr, gene){
    # EXTRACT DATA BY GENE ---
    d_tf <- df[which(df$Gene == gene),]
    d_expr <- dat_expr[which(dat_expr$Gene == gene),]

    # MERGE DATA ---
    dm <- merge(x=d_tf, y=d_expr, by=c("SampleID","Gene"), all=TRUE)

    # KEY ---
    key <- c("0:0","1:0","0:1","1:1")

    # GET MISSING KEY ---
    key_0 <- setdiff(key, unique(dm$key))
    
    # IMPUTE MISSING DATA ---
    if(length(key_0) != 0){
        dm_0 <- data.frame(SampleID="X", Gene=gene,
                            Status_znf=NA, Status_myc=NA, 
                            key=key_0, Expr=NA)
    
        dm <- rbind(dm, dm_0)
    }

    # FACTORIZE ---
    dm$key <- factor(dm$key, levels=key)

    return(dm)
}

### FUNCTION: getStats() ---
getStats <- function(dm, gene, key0, key1){
    # CHARACTER ---
    dm$key <- as.character(dm$key)

    # GET EXPR VALUES BY KEY ---
    y_0 <- dm$Expr[which(dm$key == key0)]
    y_1 <- dm$Expr[which(dm$key == key1)]

    # COMPUTE MEDIAN ---
    median_0 <- median(y_0)
	median_1 <- median(y_1)

    # COMPUTE STANDARD DEVIATION ---
    sd_0 <- sd(y_0)
	sd_1 <- sd(y_1)

    if( !is.na(median_0) & !is.na(median_1) ){
        if( !is.na(sd_0) & !is.na(sd_1) ){
            if( (sd_0 != 0) & (sd_1 != 0) ){
                if( (median_0 != 0) | (median_1 != 0) ){
                    # COMPUTE FOLDCHANGE ---
                    foldchange <- ifelse(median_1 >= median_0, median_1/median_0, -median_0/median_1)

	                # WILCONX RANK SUM TEST ---
	                res_wilcoxtest <- wilcox.test(y_0, y_1, alternative="two.sided", na.action=na.omit)

                    # GET PVALUE ---
                    pvalue <- res_wilcoxtest$p.value

	                if(pvalue == 0){
		                pvalue <- .Machine$double.eps #smallest value
	                }

                    # PREPARE DATA ---
                    d <- data.frame(Gene=gene, group0=key0, group1=key1, foldchange=foldchange, pvalue=pvalue)
                }else{
                    d <- data.frame(Gene=gene, group0=key0, group1=key1, foldchange=NA, pvalue=NA) 
                }
            }else{
                d <- data.frame(Gene=gene, group0=key0, group1=key1, foldchange=NA, pvalue=NA)   
            }
        }else{
            d <- data.frame(Gene=gene, group0=key0, group1=key1, foldchange=NA, pvalue=NA)   
        }
    }else{
        d <- data.frame(Gene=gene, group0=key0, group1=key1, foldchange=NA, pvalue=NA)            
    }

    return(d)
}

### FUNCTION: aggregateStats() ---
aggregateStats <- function(list.stat){
    ds <- do.call(rbind.data.frame, list.stat)
    ds <- subset(ds, !is.na(ds$foldchange))
    ds <- ds[order(ds$pvalue, decreasing=FALSE),]
    ds$fdr <- p.adjust(ds$pvalue, method = "BH", n = length(ds$pvalue))
    rownames(ds) <- NULL

    return(ds)
}


### FUNCTION: getVolcanoPlot() ---
getVolcanoPlot <- function(df=NULL, logfc.threshold, pvalue.threshold, filterby.pvalueadj, label_text=FALSE){
	yval <- -log10(pvalue.threshold)
    xval <- logfc.threshold
	cbPalette1 <- c("#bdbdbd","#e31a1c")

    # PVALUE FILTER ---
    if(filterby.pvalueadj == TRUE){
        df$nlogp <- -log10(df$fdr)
    }else{
        df$nlogp <- -log10(df$pvalue)
    }

    # GROUP GENES ---
    df$Group <- "STABLE"
    if(filterby.pvalueadj == TRUE){    
        df$Group[which((df$foldchange >= logfc.threshold) & (df$fdr <= pvalue.threshold))] <- "VARIABLE"
        df$Group[which((df$foldchange <= -logfc.threshold) & (df$fdr <= pvalue.threshold))] <- "VARIABLE"
    }else{
        df$Group[which((df$foldchange >= logfc.threshold) & (df$pvalue <= pvalue.threshold))] <- "VARIABLE"
        df$Group[which((df$foldchange <= -logfc.threshold) & (df$pvalue <= pvalue.threshold))] <- "VARIABLE"
    }

    # ADD GENE NAMES TO LABEL ---
    df$Gene_Label <- NA
    df$Gene_Label[which(df$Group == "VARIABLE")] <- df$Gene[which(df$Group == "VARIABLE")]

    # ORDER DATA ---
    df <- df[order(df$Group, decreasing=FALSE),]

    # FACTORIZE ---
    df$Group <- factor(df$Group, levels=c("STABLE","VARIABLE"))

    # PLOT ---
	p <- ggplot(df, aes(x=foldchange, y=nlogp, label=Gene_Label)) + 
			geom_hline(yintercept = yval, color="gray70", size=0.2, linetype=4, alpha=0.7) +
			geom_vline(xintercept = xval, color="gray70", size=0.2, linetype=4, alpha=0.7) +		
			geom_vline(xintercept = -xval, color="gray70", size=0.2, linetype=4, alpha=0.7) +
			geom_point(aes(fill=Group, color=Group, size=Group), shape = 21, stroke=0.2, alpha=1) +
			scale_fill_manual(values=cbPalette1) +
			scale_color_manual(values=cbPalette1) +
            scale_size_manual(values=c(0.2,0.8)) +
            scale_x_continuous(breaks=round(seq(-1.4, 1.4, by=0.2), 2) ) +
            scale_y_continuous(breaks=seq(0, 6, 1)) +
			coord_cartesian(xlim=c(-1.4, 1.4), ylim=c(0, 6)) +
			theme(
				axis.text = element_text(size = 5, color="#000000"),
				axis.title = element_text(size = 5, color="#000000"),
				strip.text = element_text(size = 5, color="#000000"),
				strip.background = element_rect(fill = "#FFFFFF", color = "#FFFFFF"),
				plot.title = element_text(size = 5, color="#000000", hjust=0.5),
				panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				axis.ticks = element_line(size=0.2, color="#000000"),	
				panel.background = element_rect(fill="white", color="#000000"),
				legend.text = element_text(size = 5, color="#000000"),
			    legend.title = element_blank(),
			    legend.key.size = unit(0.3, "cm"),			
			    legend.position="none") +
			xlab("log2 Fold Change") + 
			ylab("-log10(pvalue)") + 
			ggtitle("") 	

    if(label_text){
        p <- p + geom_text(aes(label=Gene_Label, alpha=0.3, hjust=1, vjust=-1, lineheight=1.5), color="#000000", size=0.5)
    }

    return(p)
}


############################################################################################################################################
### LOAD ATAC-SEQ MASTER DATA ---
list.rds_atacseq_masterdata <- readRDS(file=file.rds_atacseq_masterdata)

### GET METADATA ---
metadata <- list.rds_atacseq_masterdata$metadata_mcrpc

### LOAD GENELIST ---
list.rds_genelist <- readRDS(file=file.rds_genelist)
genes_hallmark_myc <- unique( c( list.rds_genelist$msigdb_pathways$hallmark$HALLMARK_MYC_TARGETS_V1, list.rds_genelist$msigdb_pathways$hallmark$HALLMARK_MYC_TARGETS_V2) )



################################################################################################################################################
### LOAD ANNOTATION DATA ---
annot <- parseGTF(file.annotation_gtf, feature.type="transcript")

### LOAD GENE ANNOTATION DATABASE FROM GFF/GTF ---
genecode.txdb <- GenomicFeatures::makeTxDbFromGFF(file=file.annotation_gtf, format="gtf", organism="Homo sapiens", chrominfo=seqinfo(TxDb.Hsapiens.UCSC.hg38.knownGene))



################################################################################################################################################
### LOAD FOOTPRINT DATA: BOUND AND UNBOUND ---
list.footprints_bound <- readRDS(file.rds_footprints_bound)
list.footprints_unbound <- readRDS(file.rds_footprints_unbound)


### PREPARE GENE MATRIX DATA ---
mat_genes_znf <- prepareMatrixData(list.footprints_bound, list.footprints_unbound, metadata, annot, genecode.txdb, motif_prefix=motif_prefix_znf)
mat_genes_myc <- prepareMatrixData(list.footprints_bound, list.footprints_unbound, metadata, annot, genecode.txdb, motif_prefix=motif_prefix_myc)


### EXTRACT MATRIX GENES BY GENELIST ----
mat_genes_znf_genelist <- subset(mat_genes_znf, rownames(mat_genes_znf) %in% genes_hallmark_myc)
mat_genes_myc_genelist <- subset(mat_genes_myc, rownames(mat_genes_myc) %in% genes_hallmark_myc)



################################################################################################################################################
### RESHAPE DATA ---
df_znf <- reshape2::melt(t(mat_genes_znf_genelist))
colnames(df_znf) <- c("SampleID","Gene","Status")

df_myc <- reshape2::melt(t(mat_genes_myc_genelist))
colnames(df_myc) <- c("SampleID","Gene","Status")

### MERGE DATA ---
df <- merge(x=df_znf, y=df_myc, by=c("SampleID","Gene"), all=TRUE)
colnames(df) <- c("SampleID","Gene","Status_znf","Status_myc")

### CONVERT NAs TO 0 ----------
df$Status_znf[which(is.na(df$Status_znf))] <- 0
df$Status_myc[which(is.na(df$Status_myc))] <- 0

### ORDER DATA BY GENE ---
df <- df[order(df$Gene, df$SampleID),]
rownames(df) <- NULL

### ADD KEY ---
df$key <- apply(df, 1, function(x) paste(x[3],x[4], sep=":") )

### FIX CHARACTER ---
df$SampleID <- as.character(df$SampleID)
df$Gene <- as.character(df$Gene)


### GET KEY COUNTS  ---
dm_stat <- df %>%
            dplyr::select(Gene, key) %>%
            dplyr::group_by(Gene) %>%
            dplyr::count(key, sort=FALSE, name="Freq")

dm_stat$Gene <- as.character(dm_stat$Gene)

### MATRIX BY GENE-KEYS ---
keys <- c("0:0","1:0","0:1","1:1")
mat_stat <- matrix(0, nrow=length(unique(dm_stat$Gene)), ncol=length(keys), dimnames=list(sort(unique(dm_stat$Gene)), keys))

for(i in 1:nrow(dm_stat)){
    x <- dm_stat$Gene[i]
    y <- dm_stat$key[i]
    mat_stat[x,y] <- dm_stat$Freq[i]
}

### FILTER MATRIX ----
mat_stat_pass <- mat_stat[which(mat_stat[,"1:0"] >= 3),]



################################################################################################################################################
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

### RESHAPE EXPRESSION DATA ---
dat_expr <- reshape2::melt(t(expr))
colnames(dat_expr) <- c("SampleID","Gene","Expr")

### FIX CHARACTER ---
dat_expr$SampleID <- as.character(dat_expr$SampleID)
dat_expr$Gene <- as.character(dat_expr$Gene)
######################################################################




################################################################################################################################################
### EXTRACT EXPRESSION ---
list.dm <- list()
for(i in 1:nrow(mat_stat_pass)){
    gene <- rownames(mat_stat_pass)[i]
    list.dm[[gene]] <- prepareDataStats(df, dat_expr, gene)
}


### GET STAT ---
list.stat0 <- list()
for(i in 1:length(list.dm)){
    gene <- names(list.dm)[i]
    dm <- list.dm[[gene]]
    list.stat0[[gene]] <- getStats(dm, gene, key0="0:0", key1="0:1")
}

list.stat1 <- list()
for(i in 1:length(list.dm)){
    gene <- names(list.dm)[i]
    dm <- list.dm[[gene]]
    list.stat1[[gene]] <- getStats(dm, gene, key0="0:0", key1="1:0")
}

list.stat2 <- list()
for(i in 1:length(list.dm)){
    gene <- names(list.dm)[i]
    dm <- list.dm[[gene]]
    list.stat2[[gene]] <- getStats(dm, gene, key0="0:0", key1="1:1")
}

### AGGREGATE DATA ---
dm_stat0 <- aggregateStats(list.stat=list.stat0)
dm_stat1 <- aggregateStats(list.stat=list.stat1)
dm_stat2 <- aggregateStats(list.stat=list.stat2)


################################################################################################################################################
### VOLCANO PLOT ---
p1 <- getVolcanoPlot(df=dm_stat1, logfc.threshold=1, pvalue.threshold=0.05, filterby.pvalueadj=FALSE, label_text=FALSE)
p2 <- getVolcanoPlot(df=dm_stat1, logfc.threshold=1, pvalue.threshold=0.05, filterby.pvalueadj=FALSE, label_text=TRUE)

### PLOT ---
file.plot <- file.path(dir.reproduce_fig, "figure_6_e.pdf")
pdf(file.plot, width=2, height=2)
    grid.arrange(p1, nrow=1, ncol=1)
    grid.arrange(p2, nrow=1, ncol=1)    
dev.off()









################################################################################################################################################
################ FIGURE 6F #################################################
################ BOX-PLOT EXPR (0:1, 1:0, 1:1) #############################
dm_stat_3box <- rbind(dm_stat0, dm_stat1, dm_stat2)
dm_stat_3box <- subset(dm_stat_3box, dm_stat_3box$foldchange > 0)
dm_stat_3box$group1 <- factor(dm_stat_3box$group1, levels=c("0:1","1:0","1:1"))

### PLOT ---
p_box3 <- ggplot(dm_stat_3box, aes(x=group1, y=foldchange)) +  
            geom_jitter(width=0.1, cex=0.2, color="#000000", alpha=0.5, na.rm=TRUE) +  
            geom_boxplot(aes(fill=group1), color="#000000", alpha=0.6, lwd=0.2, outlier.size=0.1, outlier.color="#bdbdbd", na.rm=FALSE) +
            scale_fill_manual(values=c("#BDBDBD","#BDBDBD","#BDBDBD")) +        
		    theme(
			    axis.text.x = element_text(size = 5, color="#000000"),
			    axis.text.y = element_text(size = 5, color="#000000"),
			    axis.title = element_text(size = 5, color="#000000"),
			    plot.title = element_text(size = 5, color="#000000", hjust=0.5),
			    panel.grid.major = element_blank(),
			    panel.grid.minor = element_blank(),
			    axis.ticks = element_line(size = 0.2, color="#000000"),	
                strip.text.x = element_text(size=5, color="#000000"),
                strip.text.y = element_text(size=5, color="#000000", angle=0, hjust=0),
                strip.background = element_rect(fill="#FFFFFF", color="#FFFFFF"),
                panel.background = element_rect(fill="#FFFFFF", color="#000000"),
			    legend.text = element_text(size = 5, color="#000000"),
			    legend.title = element_blank(),
			    legend.key.size = unit(0.2, "cm"),			
			    legend.position="none",
                legend.box = "horizontal") + 
            guides(fill = guide_legend(nrow = 2)) +                   
		    ylab("foldchange") +            
		    xlab("") + 
            ggtitle("") 

### PLOT ---
file.plot <- file.path(dir.reproduce_fig, "figure_6_f.pdf")
pdf(file.plot, width=1.5, height=2)
    grid.arrange(p_box3, nrow=1, ncol=1)  
dev.off()

### WILCOXN TEST ----
wilcox.test(dm_stat_3box$foldchange[which(dm_stat_3box$group1 == "0:1")], dm_stat_3box$foldchange[which(dm_stat_3box$group1 == "1:0")], alternative="two.sided", na.action=na.omit, paired=FALSE)$p.value
#[1] 0.313606

wilcox.test(dm_stat_3box$foldchange[which(dm_stat_3box$group1 == "0:1")], dm_stat_3box$foldchange[which(dm_stat_3box$group1 == "1:1")], alternative="two.sided", na.action=na.omit, paired=FALSE)$p.value
#[1] 0.01936213

wilcox.test(dm_stat_3box$foldchange[which(dm_stat_3box$group1 == "1:0")], dm_stat_3box$foldchange[which(dm_stat_3box$group1 == "1:1")], alternative="two.sided", na.action=na.omit, paired=FALSE)$p.value
#[1] 0.0033167

