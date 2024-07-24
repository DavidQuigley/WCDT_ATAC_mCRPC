###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### LOAD LIBRARIES ---
library("stringr")
library("data.table")
library("GenomicRanges")
library("ChIPseeker")
library("dplyr")
library("reshape2")
library("ggplot2")
library("gridExtra")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")

### DEFINE PATH ---
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

dir.footprints <- file.path(dir.reproduce_data, "tf_footprints/mcrpc_subtypes")


### DEFINE FILES ---
file.rds_footprints <- file.path(dir.footprints, "data_footprints/tobias_mcrpc_subtypes_tf_footprints_featureids_bound_unbound.rds")
file.rds_tf_hits <- file.path(dir.footprints, "data_footprints/tobias_mcrpc_subtypes_tf_hits.rds")
file.annotation_gtf <- file.path(dir.wrk, "reference/annotation_data/gene_annotation/gencode_v28/gencode.v28.annotation.gtf.gz")


### FUNCTION: getGRanges() ---
getGRanges <- function(feature_ids){
    # GET BED FILE ---
    df <- data.frame(chr=unlist(lapply(stringr::str_split(feature_ids, "_"), function(x) x[1])),
                    start=as.numeric(unlist(lapply(stringr::str_split(feature_ids, "_"), function(x) x[2]))),
                    end=as.numeric(unlist(lapply(stringr::str_split(feature_ids, "_"), function(x) x[3]))))

    # GET GRANGES ---
    gr <- GenomicRanges::makeGRangesFromDataFrame(df, keep.extra.columns=FALSE)
    gr <- sort(gr)
    gr <- gr[unique(findOverlaps(gr, type = "any", select = "first"))]

    return(gr)
}

### FUNCTION: annotPeaksbyRange() ---
annotPeaksbyRange <- function(gr_peaks, genecode.txdb, motif_prefix){
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

    # ADD SAMPLEID --
    annot_peaks$motif_prefix <- motif_prefix

    return(annot_peaks)
}

### LOAD FOOTPRINTS ---
list.rds_footprints <- readRDS(file=file.rds_footprints)

### LOAD TF HITS ---
list.rds_tf_hits <- readRDS(file=file.rds_tf_hits)
tfs <- list.rds_tf_hits$tf_footprints_hits_score$motif_prefix

### LOAD GENE ANNOTATION DATABASE FROM GFF/GTF ---
genecode.txdb <- GenomicFeatures::makeTxDbFromGFF(file=file.annotation_gtf, format="gtf", organism="Homo sapiens", chrominfo=seqinfo(TxDb.Hsapiens.UCSC.hg38.knownGene))

### GET BOUND FOOTPRINTS GRANGES ---
list.gr_motif <- list()
for(motif_prefix in tfs){
    # GET MOTIF FEATURES ---
    feature_ids <- unique( unlist(list.rds_footprints$bound[[motif_prefix]]) )

    # GET MOTIF RANGES ---
    list.gr_motif[[motif_prefix]] <- getGRanges(feature_ids)
}


### GET MOTIF POSITION ANNOTATION ---
list.annot <- list()
for(motif_prefix in tfs){
    list.annot[[motif_prefix]] <- annotPeaksbyRange(gr_peaks=list.gr_motif[[motif_prefix]], genecode.txdb, motif_prefix)
    cat("PROCESSED", motif_prefix, "\n", sep=" ")
}




###########################################################################################################################
### EXTRACT COLUMNS ---
list.dat <- lapply(list.annot, function(x) subset(x, select=c("motif_prefix","geneId","Feature","distanceToTSS")))

### AGGREGATE DATA ---
dat <- do.call(rbind.data.frame, list.dat)
rownames(dat) <- NULL

### FORMAT DIST TO TSS --
dat$distanceToTSS_kb <- abs(dat$distanceToTSS/1000)
dat$distanceToTSS_kb_log <- log2(dat$distanceToTSS_kb + 1)

### GET STATS ---
df <- dat %>% 
        dplyr::select(motif_prefix, distanceToTSS_kb_log) %>%
        dplyr::group_by(motif_prefix) %>%
        dplyr::summarize(Quartile1=round( as.numeric(summary(distanceToTSS_kb_log)[2]) , 3),
                        Mean=round(mean(distanceToTSS_kb_log), 3), 
                        Median=round(median(distanceToTSS_kb_log), 3), 
                        SD=round(sd(distanceToTSS_kb_log), 3)) %>%
        dplyr::arrange(Mean)

### SEPARATE DATA ---
df1 <- df[102:203,]
df2 <- df[1:101,]

### PREPARE DATA ---
dat1 <- dat[which(dat$motif_prefix %in% df1$motif_prefix),]
dat2 <- dat[which(dat$motif_prefix %in% df2$motif_prefix),]

### FACTORIZE ---
dat1$motif_prefix <- factor(dat1$motif_prefix, levels=df1$motif_prefix)
dat2$motif_prefix <- factor(dat2$motif_prefix, levels=df2$motif_prefix)



### FUNCTION: getBoxPlot() ---
getBoxPlot <- function(dm){
    # PLOT ---
    p <- ggplot(dm, aes(x=distanceToTSS_kb_log, y=motif_prefix)) +
            geom_boxplot(fill="#fdbf6f", color="#000000", alpha=0.6, lwd=0.2, outlier.size=0.1, outlier.color="#bdbdbd", na.rm=TRUE) +
            coord_cartesian(xlim=c(0,12)) +
            scale_x_continuous(breaks=seq(0,12, by=2), labels=c(0,3,15,63,255,1023,4095)) +
            geom_vline(xintercept=2, color="#e41a1c", size=0.8, linetype=2) +
		    theme(
			    axis.text.x = element_text(size = 4.5, color="#000000"),
			    axis.text.y = element_text(size = 4.5, color="#000000"),
			    axis.title = element_text(size = 5, color="#000000"),
			    plot.title = element_text(size = 5, color="#000000", hjust=0.5),
			    panel.grid.major = element_blank(),
			    panel.grid.minor = element_blank(),
			    axis.ticks = element_line(size = 0.2, color="#000000"),	
			    panel.background = element_rect(fill = "#FFFFFF", color = "#000000"),
			    legend.text = element_text(size = 5, color="#000000"),
			    legend.title = element_blank(),
			    legend.key.size = unit(0.2, "cm"),			
			    legend.position="none",
                legend.box = "horizontal") + 
            guides(fill = guide_legend(nrow = 2)) +                   
		    ylab("") +            
		    xlab("absolute distance to TSS in kb") + 
            ggtitle("") 

    return(p)
}


# PLOT ---
p1 <- getBoxPlot(dm=dat1)
p2 <- getBoxPlot(dm=dat2)

### GENERATE PLOT ---
file.plot <- file.path(dir.reproduce_fig, "supplementary_figure_17.pdf")
pdf(file.plot, width=6.8, height=8)
    grid.arrange(p1, p2, nrow=1, ncol=2)
dev.off()


