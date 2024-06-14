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

### FUNCTION: computeStats() ---
computeStats <- function(dm, gene){
    # GET VALUES ---
    y_bound <- dm$Expr[which(dm$Group == "BOUND")]
    y_unbound <- dm$Expr[which(dm$Group == "UNBOUND")]

    # COMPUTE MEDIAN ---
    median_bound <- median(y_bound)
	median_unbound <- median(y_unbound)

    # COMPUTE STANDARD DEVIATION ---
    sd_bound <- sd(y_bound)
	sd_unbound <- sd(y_unbound)

    if( (sd_bound != 0) & (sd_unbound != 0) ){
        if( (median_bound != 0) | (median_unbound != 0) ){
            # COMPUTE FOLDCHANGE ---
            foldchange <- ifelse(median_bound >= median_unbound, median_bound/median_unbound, -median_unbound/median_bound)

	        # WILCONX RANK SUM TEST ---
	        res_wilcoxtest <- wilcox.test(y_bound, y_unbound, alternative="two.sided", na.action=na.omit)

            # GET PVALUE ---
            pvalue <- res_wilcoxtest$p.value

	        if(pvalue == 0){
		        pvalue <- .Machine$double.eps #smallest value
	        }

            # PREPARE DATA ---
            d <- data.frame(Gene=gene, foldchange=foldchange, pvalue=pvalue)
        }else{
            d <- data.frame()            
        }
    }else{
        d <- data.frame()
    }

    return(d)
}



### FUNCTION: getVolcanoPlot() ---
getVolcanoPlot <- function(file=NULL, df=NULL, logfc.threshold, pvalue.threshold, filterby.pvalueadj, fc.upperbound, nlogp.upperbound, label_text=FALSE){
	yval <- -log10(pvalue.threshold)
    xval <- logfc.threshold
	cbPalette1 <- c("#bdbdbd","#1f78b4","#e31a1c")

    # LOAD DEA SUMMARY ---
    if(!is.null(file)){
        df <- data.table::fread(file=file, sep="\t", header=TRUE, nThread=50, data.table=FALSE, verbose=FALSE)
    }

    # SUBSET DATA BY COLUMN ---
    #df <- subset(df, select=c("FeatureID","log2FoldChange","pvalue","padj","Gene"))
    
    # REMOVE NAs ---
    #if(length(which(is.na(df$pvalue))) > 0){
    #    df <- df[-which(is.na(df$pvalue)),]
    #}

    # PVALUE FILTER ---
    if(filterby.pvalueadj == TRUE){
        df$nlogp <- -log10(df$fdr)
    }else{
        df$nlogp <- -log10(df$pvalue)
    }

    # GROUP GENES ---
    df$Group <- "STABLE"
    if(filterby.pvalueadj == TRUE){    
        df$Group[which((df$foldchange >= logfc.threshold) & (df$fdr <= pvalue.threshold))] <- "VARIABLE_POS"
        df$Group[which((df$foldchange <= -logfc.threshold) & (df$fdr <= pvalue.threshold))] <- "VARIABLE_NEG"
    }else{
        df$Group[which((df$foldchange >= logfc.threshold) & (df$pvalue <= pvalue.threshold))] <- "VARIABLE_POS"
        df$Group[which((df$foldchange <= -logfc.threshold) & (df$pvalue <= pvalue.threshold))] <- "VARIABLE_NEG"
    }

    # LABEL GROUP BY PATHWAY GENES TO HIGHLIGHT ---
    #df$Group[which((df$Group == "VARIABLE") & (df$Gene %in% genelist))] <- "ZHIGHLIGHT"

    # ADD GENE NAMES TO LABEL ---
    df$Gene_Label <- NA
    df$Gene_Label[which(df$Group %in% c("VARIABLE_POS","VARIABLE_NEG"))] <- df$Gene[which(df$Group %in% c("VARIABLE_POS","VARIABLE_NEG"))]

    # ADJUST FOLDCHANGE FOR VIZ ---
    if(max(df$foldchange) > fc.upperbound){
        df$foldchange[which(df$foldchange > fc.upperbound)] <- fc.upperbound 
    }

    if(min(df$foldchange) < -fc.upperbound){
        df$foldchange[which(df$foldchange < -fc.upperbound)] <- -fc.upperbound
    }

    # ADJUST nlogp FOR VIZ ---
    if(max(df$nlogp) > nlogp.upperbound){
        df$nlogp[which(df$nlogp > nlogp.upperbound)] <- nlogp.upperbound 
    }

    # ORDER DATA ---
    df <- df[order(df$Group, decreasing=FALSE),]

    # FACTORIZE ---
    df$Group <- factor(df$Group, levels=c("STABLE","VARIABLE_NEG","VARIABLE_POS"))

    # PLOT ---
	p <- ggplot(df, aes(x=foldchange, y=nlogp, label=Gene_Label)) + 
			geom_hline(yintercept = yval, color="gray70", size=0.2, linetype=4, alpha=0.7) +
			geom_vline(xintercept = xval, color="gray70", size=0.2, linetype=4, alpha=0.7) +		
			geom_vline(xintercept = -xval, color="gray70", size=0.2, linetype=4, alpha=0.7) +
			geom_point(aes(fill=Group, color=Group, size=Group), shape = 21, stroke=0.2, alpha=1) +
			scale_fill_manual(values=cbPalette1) +
			scale_color_manual(values=cbPalette1) +
            scale_size_manual(values=c(0.2,0.8,0.8)) +
            scale_x_continuous(breaks=seq(-2, fc.upperbound, 1)) +
            scale_y_continuous(breaks=seq(0, 4, 1)) +
			coord_cartesian(xlim=c(-2, fc.upperbound), ylim=c(0, 4)) +
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


################################################################################################################################################
### LOAD ANNOTATION DATA ---
annot <- parseGTF(file.annotation_gtf, feature.type="transcript")

### LOAD GENE ANNOTATION DATABASE FROM GFF/GTF ---
genecode.txdb <- GenomicFeatures::makeTxDbFromGFF(file=file.annotation_gtf, format="gtf", organism="Homo sapiens", chrominfo=seqinfo(TxDb.Hsapiens.UCSC.hg38.knownGene))




################################################################################################################################################
### LOAD FOOTPRINT DATA: BOUND ---
list.footprints_bound <- readRDS(file.rds_footprints_bound)
list.footprints_bound <- list.footprints_bound[[motif_prefix_znf]]
names(list.footprints_bound) <- metadata$Sample_ID

### LOAD FOOTPRINT DATA: UNBOUND ---
list.footprints_unbound <- readRDS(file.rds_footprints_unbound)
list.footprints_unbound <- list.footprints_unbound[[motif_prefix_znf]]
names(list.footprints_unbound) <- metadata$Sample_ID

### GET FEATUREIDS ---
feature_ids_bound <- unique(unlist(lapply(list.footprints_bound, function(x) x$motif_position)))
feature_ids_unbound <- unique(unlist(lapply(list.footprints_unbound, function(x) x$motif_position)))

### GET UNIQUE FEATUREIDS ---
feature_ids <- unique(c(feature_ids_bound, feature_ids_unbound))
gr_footprints <- getFootprintRanges(feature_ids)


### GET FOOTPRINT SITE ANNOTATION -----
dat_footprint_annot <- annotPeaksbyRange(gr_peaks=gr_footprints, genecode.txdb)
dat_footprint_annot <- subset(dat_footprint_annot, dat_footprint_annot$Feature == "Promoter")

### SUBSET FEATUREIDS MATRIX WITH FOOTPRINTS ON PROMOTERS ---
dat_footprint_annot <- addGenes(annot, dat=dat_footprint_annot)
dat_footprint_annot_prom <- subset(dat_footprint_annot, dat_footprint_annot$GeneType == "protein_coding")


### GET GRANGES ---
list.footprints_bound_gr <- lapply(list.footprints_bound, function(x) GenomicRanges::makeGRangesFromDataFrame(x, keep.extra.columns=TRUE) )
list.footprints_unbound_gr <- lapply(list.footprints_unbound, function(x) GenomicRanges::makeGRangesFromDataFrame(x, keep.extra.columns=TRUE) )



################################################################################################################################################
### CREATE GENE LEVEL MATRIX FOR BOUND/UNBOUND IN PROMOTERS ---
genes <- sort(unique(dat_footprint_annot_prom$Gene))
mat_genes <- matrix(0, nrow=length(genes), ncol=length(list.footprints_bound_gr), dimnames=list(genes, names(list.footprints_bound_gr)))

### LABEL BOUND GENES IN THE MATRIX ----
for(i in 1:length(list.footprints_bound_gr)){
    sampleid <- names(list.footprints_bound_gr)[i]
    gr <- list.footprints_bound_gr[[sampleid]]
    f_ids <- gr$motif_position

    genes_bound <- unique(dat_footprint_annot_prom$Gene[which(dat_footprint_annot_prom$FeatureID %in% f_ids)])

    index_bound <- which(rownames(mat_genes) %in% genes_bound)
    mat_genes[index_bound,i] <- 1

    cat("PROCESSED:", sampleid, "\n", sep=" ")
}

#> dim(mat_genes)
#[1] 16509    70

### REMOVE GENES TF-UNBOUND/BOUND IN ALL SAMPLES OR JUST 1 SAMPLE ----
mat_genes <- mat_genes[-which(rowSums(mat_genes) == 0),]
mat_genes <- mat_genes[-which(rowSums(mat_genes) == length(list.footprints_bound_gr)),]
mat_genes <- mat_genes[-which(rowSums(mat_genes) == 1),]
mat_genes <- mat_genes[-which(rowSums(mat_genes) == (length(list.footprints_bound_gr)-1) ),]

#> dim(mat_genes)
#[1] 6331   70


### MELT DATA ---
df_genes <- reshape2::melt(t(mat_genes))
colnames(df_genes) <- c("SampleID","Gene","Status")
df_genes$Group <- ifelse(df_genes$Status , "BOUND", "UNBOUND")





############################################################################################################################################
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
######################################################################

### SUBSET BY EXPRESSION DATA GENES ---
expr <- subset(expr, rownames(expr) %in% rownames(mat_genes))
expr <- expr[match(rownames(mat_genes), rownames(expr)),]

### MELT DATA ---
df_expr <- reshape2::melt(t(expr))
colnames(df_expr) <- c("SampleID","Gene","Expr")

### ADD EXPRESSION DATA ---
df_genes$Expr <- df_expr$Expr


############################################################################################################################################

########## RESTRICT ANALYSIS BY GENE EXPRESSION LEVELS -----------------
### INCLUDE GENES WITH AT LEAST ONE SAMPLE EXPRESSING LOG2 TPM >= 5
genes <- as.character(unique(df_genes$Gene))
list.genes_high <- list()
ctr <- 1
for(i in 1:length(genes)){
	gene <- genes[i]
	dm <- subset(df_genes, df_genes$Gene == gene)

    if(length(which(dm$Expr >= 5)) != 0){
        list.genes_high[[ctr]] <- gene
        ctr <- ctr + 1
    }

    cat("PROCESSED:", i, "\n", sep=" ")
}

### GENES HIGH ---
genes_high <- unlist(list.genes_high)




### COMPUTE STATS ---
list.df <- list()
for(i in 1:length(genes_high)){
	gene <- genes_high[i]
	dm <- subset(df_genes, df_genes$Gene == gene)
    list.df[[gene]] <- computeStats(dm, gene)
    cat("PROCESSED:", i, "\n", sep=" ")
}

### AGGREGATE DATA ---
df <- do.call(rbind.data.frame, list.df)
rownames(df) <- NULL
df <- df[-which(df$foldchange %in% c("Inf","-Inf") ),]

### COMPUTE FDR ---
df <- df[order(df$pvalue, decreasing=FALSE),]
df$fdr <- p.adjust(df$pvalue, method = "BH", n = length(df$pvalue))
rownames(df) <- NULL




############################################################################################################################################

### GET VOLCANO PLOT ---
p1 <- getVolcanoPlot(file=NULL, df=df, logfc.threshold=1, pvalue.threshold=0.05, filterby.pvalueadj=TRUE, fc.upperbound=6, nlogp.upperbound=4, label_text=FALSE)
p2 <- getVolcanoPlot(file=NULL, df=df, logfc.threshold=1, pvalue.threshold=0.05, filterby.pvalueadj=TRUE, fc.upperbound=6, nlogp.upperbound=4, label_text=TRUE)


### PLOT ---
file.plot <- file.path(dir.reproduce_fig, "figure_6_a.pdf")
pdf(file.plot, width=2, height=2)
    grid.arrange(p1, nrow=1, ncol=1)
    grid.arrange(p2, nrow=1, ncol=1)    
dev.off()




############################################################################################################################################
############################################ OVER REPRESNTATION ANALYSIS ###################################################################

### EXTRACT UPREGULATED GENES UPON ZNF263 BINDING ---
df_up <- df[which( (df$foldchange > 0) & (df$fdr < 0.05 ) ),]

### DEFINE GENES ---
list.genesets <- list.rds_genelist$msigdb_pathways$hallmark
genes.queryset <- df_up$Gene
genes.refset <- gannot$Gene  #BACKGROUND GENESET: ALL GENES IN RNA-SEQ DATA 


# GET OVER-REPRESENTATION ANALYSIS (HYPER-GEOMERTIC TEST) ---
dat_ora <- get_over_representation_analysis(list.genesets, genes.queryset, genes.refset, p.threshold=0.0001)





############################################################################################################################################
############################################### PLOT ENRICHMENT ############################################################################

### FUNCTION: getPathwaysData() ---
getPathwaysData <- function(dat){
    # EXTRACT DATA ---
    if(nrow(dat) >= 15){
        dat <- dat[1:15,]
    }

    # ADD PVALUES ---
    pval <- min(setdiff(dat$pvalue, 0))
    index0 <- which(dat$pvalue == 0)

    if(length(index0) != 0){
        dat$pvalue[which(dat$pvalue == 0)] <- pval
    }

    # ADD VALUES ---
    dat$nlogp <- -log10(dat$fdr)
    dat <- dat[order(dat$nlogp, decreasing=FALSE),]
    dat$Category <- factor(dat$Category, levels=dat$Category)

    return(dat)
}

### FUNCTION: plotBar() ----
plotBar <- function(df){
    # PLOT ---
    p <- ggplot(df, aes(x=nlogp, y=Category)) +
            geom_bar(stat="identity", fill="#006094") +
            coord_cartesian(xlim=c(0,30)) +
            scale_x_continuous(breaks=seq(0,30,by=5)) +
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
            xlab("-log10(fdr)") +
            ylab("") + 
            ggtitle("") 
    
    return(p)
}


### LOAD ENRICHMENT DATA ---
dat_ora_parsed <- getPathwaysData(dat=dat_ora)

### GET PLOT ---
p_bar <- plotBar(df=dat_ora_parsed)

### PLOT ---
file.plot <- file.path(dir.reproduce_fig, "supplementary_figure_22_b.pdf")
pdf(file.plot, width=2.5, height=2)
    grid.arrange(p_bar, nrow=1, ncol=1)
dev.off()



