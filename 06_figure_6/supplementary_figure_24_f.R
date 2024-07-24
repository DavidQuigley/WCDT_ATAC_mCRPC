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
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

dir.footprints <- file.path(dir.reproduce_data, "tf_footprints/mcrpc_subtypes")

### DEFINE FILES ---
file.rds_atacseq_masterdata <- file.path(dir.reproduce_data, "wcdt_atacseq_masterdata.rds")
file.rds_rnaseq <- file.path(dir.reproduce_data, "wcdt_rnaseq_atacseq_gene_expression.rds")
file.rds_genelist <- file.path(dir.reproduce_data, "genelist.rds")
file.rds_footprints <- file.path(dir.footprints, "data_footprints/tobias_mcrpc_subtypes_tf_footprints_featureids_bound_unbound.rds")
file.annotation_gtf <- file.path(dir.wrk, "reference/annotation_data/gene_annotation/gencode_v28/gencode.v28.annotation.gtf.gz")
file.utility_functions <- file.path(dir.reproduce_scripts, "00_prepare_data/00_01_utility_functions.R")

### LOAD FUNCTIONS ---
source(file.utility_functions)

### SUBTYPES ---
subtypes <- c("ARpNEn","ARnNEn")

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


### FUNCTION: getFootprintMatrix() ---
getFootprintMatrix <- function(gr_footprints, list.footprints_bound_motif, subtypes){
    # CONVERT TO DATAFRAME ---
    d <- as.data.frame(gr_footprints)

    # GET featureids ---
    feature_ids <- d$FeatureID

    # CREATE MATRIX ---
    mat_features <- matrix(0, nrow=length(feature_ids), ncol=length(subtypes), dimnames=list(feature_ids, subtypes))

    # FILL MATRIX ---
    for(subtype in subtypes){
        fids_bound <- list.footprints_bound_motif[[subtype]]
        mat_features[which(feature_ids %in% fids_bound), subtype] <- 1
    }

    return(mat_features)
}

### FUNCTION: getGeneLevelMatrix() ---
getGeneLevelMatrix <- function(gr_footprints, list.footprints_bound_motif, dat_annot, subtypes){
    # CONVERT TO DATAFRAME ---
    d <- as.data.frame(gr_footprints)

    # GET GENES ---
    genes <- sort(unique(dat_annot$Gene))

    # CREATE MATRIX ---
    mat_genes <- matrix(0, nrow=length(genes), ncol=length(subtypes), dimnames=list(genes, subtypes))

    # FILL MATRIX ---
    for(subtype in subtypes){
        fids_bound <- list.footprints_bound_motif[[subtype]]
        genes_bound <- sort( unique( dat_annot$Gene[which(dat_annot$FeatureID %in% fids_bound)] ) )
        mat_genes[which(rownames(mat_genes) %in% genes_bound), subtype] <- 1
    }

    return(mat_genes)
}



### FUNCTION: prepareMatrixData() ---
prepareMatrixData <- function(list.footprints_bound, list.footprints_unbound, subtypes, annot, genecode.txdb, motif_prefix){
    # EXTRACT MOTIF FOOTPRINTS ----
    list.footprints_bound_motif <- list.footprints_bound[[motif_prefix]]
    list.footprints_unbound_motif <- list.footprints_unbound[[motif_prefix]]

    # EXTRACT BY SUBTYPES ---
    list.footprints_bound_motif <- list.footprints_bound_motif[subtypes]
    list.footprints_unbound_motif <- list.footprints_unbound_motif[subtypes]

    # GET FEATUREIDS ---
    feature_ids_bound <- unique(unlist(lapply(list.footprints_bound_motif, function(x) x)))
    feature_ids_unbound <- unique(unlist(lapply(list.footprints_unbound_motif, function(x) x)))

    # GET UNIQUE FEATUREIDS ---
    feature_ids <- unique(c(feature_ids_bound, feature_ids_unbound))
    gr_footprints <- getFootprintRanges(feature_ids)

    # GET FOOTPRINT SITE ANNOTATION -----
    dat_footprint_annot <- annotPeaksbyRange(gr_peaks=gr_footprints, genecode.txdb)
    dat_footprint_annot <- subset(dat_footprint_annot, dat_footprint_annot$Feature == "Promoter")

    # ADD GENE INFO ---
    dat_footprint_annot <- addGenes(annot, dat=dat_footprint_annot)
    dat_footprint_annot <- subset(dat_footprint_annot, dat_footprint_annot$GeneType == "protein_coding")

    # GET FOOTPRINT MATRIX ---
    #mat <- getFootprintMatrix(gr_footprints, list.footprints_bound_motif, subtypes)

    # GET GENE LEVEL MATRIX ---
    mat_genes <- getGeneLevelMatrix(gr_footprints, list.footprints_bound_motif, dat_annot=dat_footprint_annot, subtypes)

    return(mat_genes)
}



############################################################################################################################################
### LOAD ATAC-SEQ MASTER DATA ---
list.rds_atacseq_masterdata <- readRDS(file=file.rds_atacseq_masterdata)

### GET METADATA ---
metadata <- list.rds_atacseq_masterdata$metadata_mcrpc
metadata <- subset(metadata, metadata$mCRPC_Subtype %in% c("AR+_NE-","AR-_NE-"))
metadata <- metadata[order(metadata$mCRPC_Subtype, decreasing=TRUE),]

### SAMPLES BY SUBTYPE ---
ids_arpc <- metadata$Sample_ID[which(metadata$mCRPC_Subtype == "AR+_NE-")]
ids_dnpc <- metadata$Sample_ID[which(metadata$mCRPC_Subtype == "AR-_NE-")]

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
list.footprints <- readRDS(file.rds_footprints)
list.footprints_bound <- list.footprints$bound
list.footprints_unbound <- list.footprints$unbound


### PREPARE FOOTPRINTS MATRIX DATA ---
mat_genes_znf <- prepareMatrixData(list.footprints_bound, list.footprints_unbound, subtypes, annot, genecode.txdb, motif_prefix=motif_prefix_znf)
mat_genes_myc <- prepareMatrixData(list.footprints_bound, list.footprints_unbound, subtypes, annot, genecode.txdb, motif_prefix=motif_prefix_myc)


### EXTRACT MATRIX GENES BY GENELIST ----
mat_genes_znf_genelist <- subset(mat_genes_znf, rownames(mat_genes_znf) %in% genes_hallmark_myc)
mat_genes_myc_genelist <- subset(mat_genes_myc, rownames(mat_genes_myc) %in% genes_hallmark_myc)



################################################################################################################################################
### RESHAPE DATA ---
df_znf <- reshape2::melt(t(mat_genes_znf_genelist))
colnames(df_znf) <- c("Subtype","Gene","Status")

df_myc <- reshape2::melt(t(mat_genes_myc_genelist))
colnames(df_myc) <- c("Subtype","Gene","Status")

### MERGE DATA ---
df <- merge(x=df_znf, y=df_myc, by=c("Subtype","Gene"), all=TRUE)
colnames(df) <- c("Subtype","Gene","Status_znf","Status_myc")

### CONVERT NAs TO 0 ----------
df$Status_znf[which(is.na(df$Status_znf))] <- 0
df$Status_myc[which(is.na(df$Status_myc))] <- 0

### ORDER DATA BY GENE ---
df <- df[order(df$Gene, df$Subtype),]
rownames(df) <- NULL

### ADD KEY ---
df$key <- apply(df, 1, function(x) paste(x[3],x[4], sep=":") )

### FIX CHARACTER ---
df$Subtype <- as.character(df$Subtype)
df$Gene <- as.character(df$Gene)




################################################################################################################################################
######### BAR PLOT FOR MYC TARGETS KEY ######################
dm_key_freq <- df %>%
                dplyr::select(Subtype, key) %>%
                dplyr::group_by(Subtype) %>%
                dplyr::count(key, sort=FALSE, name="Freq")

### FACTORIZE ---
dm_key_freq$Subtype <- factor(dm_key_freq$Subtype, levels=c("ARnNEn","ARpNEn"))
dm_key_freq$key <- factor(dm_key_freq$key, levels=rev(c("0:0","1:0","0:1","1:1")))


### PLOT BAR ----
p  <- ggplot(dm_key_freq, aes(x=Freq, y=key)) +
        geom_bar(stat="identity", fill="#006094", color="#000000", width=0.8, size=0.2) +
        #coord_cartesian(ylim=c(40,100) ) +
        #scale_y_continuous(breaks=seq(40,100, by=10) ) +
        facet_wrap(~Subtype, nrow=1, ncol=2, scales="fixed") +  
        theme(
            axis.text.x = element_text(size = 5, color="#000000"),
            axis.text.y = element_text(size = 5, color="#000000"),
            axis.title = element_text(size = 5, color="#000000"),
            plot.title = element_text(size = 5, color="#000000", hjust=0.5),
            panel.grid.major.y = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor = element_blank(),
            axis.ticks.y = element_line(size=0.2, color="#000000"),
            axis.ticks.x = element_line(size=0.2, color="#000000"),
            strip.text = element_text(size=5, color="#000000"),
            strip.background = element_rect(fill="#FFFFFF", color="#FFFFFF"),
            panel.background = element_rect(fill="#FFFFFF", color="#000000"),
            legend.text = element_text(size = 5, color="#000000"),
            legend.title = element_blank(),
            legend.key.size = unit(0.5, "cm"),
            legend.position = "none") +
        xlab("No. of genes") +
        ylab("") + 
        ggtitle("HALLMARK_MYC_TARGETS_V1/V2")


### BOX PLOT EXPR ---
file.plot <- file.path(dir.reproduce_fig, "supplementary_figure_24_f.pdf")
pdf(file.plot, width=2, height=1.5)
    grid.arrange(p, nrow=1, ncol=1)
dev.off()


