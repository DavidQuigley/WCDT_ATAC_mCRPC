###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### LOAD LIBRARIES ---
library("stringr")
library("data.table")
library("GenomicRanges")
library("ChIPseeker")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("ggplot2")
library("gridExtra")

### DEFINE PATH ---
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/analysis/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

dir.footprints <- file.path(dir.reproduce_data, "tf_footprints/mcrpc_subtypes/data_footprints") 
dir.chipatlas <- file.path(dir.wrk, "external_data/chipatlas/data") # PATH CONTAINING CHIPATLAS FILES

### DEFINE FILES ---
file.annotation_gtf <- file.path(dir.wrk, "reference/annotation_data/gene_annotation/gencode_v28/gencode.v28.annotation.gtf.gz")
file.tf_hits_rds <- file.path(dir.footprints, "tobias_mcrpc_subtypes_tf_hits.rds")

### CHROMOSOMES ---
chromosomes <- paste("chr", c(1:22,"X","Y"), sep="")

### SUBTYPES --
subtypes <- c("ARpNEn","ARlNEn","ARpNEp","ARnNEp","ARnNEn")

###########################################################################################################################
### FUNCTION: get_motif_by_subtype() ---
get_motif_by_subtype <- function(dat.tfsubtypes, subtype){
    # EXTRACT TF BY SUBTYPE ---
    index_subtype <- which( colnames(dat.tfsubtypes) == subtype )
    index_tf <- which( dat.tfsubtypes[,index_subtype] >= 2 ) 
    dat.tfsubtypes <- dat.tfsubtypes[index_tf,]
    tfs_subtype <- dat.tfsubtypes$gene

    return(tfs_subtype)
}

### FUNCTION: parse_chipseq_data() ---
parse_chipseq_data <- function(file.chipseq, chromosomes){
    # LOAD DATA ---
    dat <- data.table::fread(file=file.chipseq, sep="\t", header=FALSE, nThread=30, data.table=FALSE, verbose=FALSE)
    
    # ADD COLUMN NAMES ---
    colnames(dat)[1] <- "chr"
    colnames(dat)[2] <- "start"
    colnames(dat)[3] <- "end"
    colnames(dat)[4] <- "attr"

    # EXTRACT DATA BY CHROMOSOMES ---
    dat <- subset(dat, dat$chr %in% chromosomes)

    # SPLIT ATTRIBUTES COLUMN ---
    #list.attr <- stringr::str_split(dat$attr, ";")
    #dat$ID <- unlist(lapply(stringr::str_split( unlist(lapply(list.attr, function(x) x[1])), "="), function(x) x[2]))
    #dat$Name <- unlist(lapply(stringr::str_split( unlist(lapply(list.attr, function(x) x[2])), "="), function(x) x[2]))
    #dat$GSMID <- unlist(lapply(stringr::str_split( unlist(lapply(stringr::str_split( unlist(lapply(list.attr, function(x) x[3])), "="), function(x) x[2])), ":"), function(x) x[1]))
    #dat$CellGroup <- unlist(lapply(stringr::str_split( unlist(lapply(list.attr, function(x) x[4])), "="), function(x) x[2]))
    #dat$SourceName <- unlist(lapply(stringr::str_split( unlist(lapply(list.attr, function(x) x[5])), "="), function(x) x[2]))
    #dat$FeatureID <- apply(dat, 1, function(x) paste(x[1], as.numeric(x[2]), as.numeric(x[3]), sep="_"))

    return(dat)
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

### FUNCTION: getBoxPlot() ---
getBoxPlot <- function(dm){
    # PLOT ---
    p <- ggplot(dm, aes(x=distanceToTSS_kb_log, y=motif_prefix)) +
            geom_boxplot(fill="#fdbf6f", color="#000000", alpha=0.6, lwd=0.2, outlier.size=0.1, outlier.color="#bdbdbd", na.rm=TRUE) +
            coord_cartesian(xlim=c(0,12)) +
            scale_x_continuous(breaks=seq(0,12, by=2), labels=c(0,3,15,63,255,1023,4095)) +
            geom_vline(xintercept=2, color="#e41a1c", size=0.8, linetype=2) +
		    theme(
			    axis.text.x = element_text(size = 6, color="#000000"),
			    axis.text.y = element_text(size = 6, color="#000000"),
			    axis.title = element_text(size = 8, color="#000000"),
			    plot.title = element_text(size = 8, color="#000000", hjust=0.5),
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

### FUNCTION: get_tf_peaks_distTSS() ---
get_tf_peaks_distTSS <- function(genecode.txdb, dir.chipatlas, chromosomes, tfs){
    # DEFINE FILES ---
    files.chipseq <- sapply(tfs, function(x) { paste(dir.chipatlas, "/Oth.ALL.05.", sprintf("%s", x), ".AllCell.bed", sep="") })
    files.chipseq <- files.chipseq[file.exists(files.chipseq)]
    tfs <- names(files.chipseq)

    # GET CHIPSEQ PEAKS FOR EACH TF ---
    list.chipseq <- list()
    for(tf in tfs){
        cat("START:", tf, "\n", sep=" ")
            list.chipseq[[tf]] <- parse_chipseq_data(file.chipseq=files.chipseq[tf], chromosomes=chromosomes)
        cat("PROCESSED:", tf, "\n", sep=" ")
    }

    # CONVERT TO GRANGE LIST ---
    list.gr_chipseq <- lapply(list.chipseq, function(x) GenomicRanges::makeGRangesFromDataFrame(x, keep.extra.columns=FALSE) )

    # GET MOTIF POSITION ANNOTATION ---
    list.annot <- list()
    for(motif_prefix in tfs){
        list.annot[[motif_prefix]] <- annotPeaksbyRange(gr_peaks=list.gr_chipseq[[motif_prefix]], genecode.txdb, motif_prefix)
        cat("PROCESSED", motif_prefix, "\n", sep=" ")
    }

    # AGGREGATE DATA ---
    dat <- do.call(rbind.data.frame, list.annot)
    rownames(dat) <- NULL

    # FORMAT DIST TO TSS --
    dat$distanceToTSS_kb <- abs(dat$distanceToTSS/1000)
    dat$distanceToTSS_kb_log <- log2(dat$distanceToTSS_kb + 1)

    # GET STATS ---
    df <- dat %>% 
            dplyr::select(motif_prefix, distanceToTSS_kb_log) %>%
            dplyr::group_by(motif_prefix) %>%
            dplyr::summarize(Quartile1=round( as.numeric(summary(distanceToTSS_kb_log)[2]) , 3),
                        Mean=round(mean(distanceToTSS_kb_log), 3), 
                        Median=round(median(distanceToTSS_kb_log), 3), 
                        SD=round(sd(distanceToTSS_kb_log), 3)) %>%
            dplyr::arrange(Median)

    # FACTORIZE ---
    dat$motif_prefix <- factor(dat$motif_prefix, levels=df$motif_prefix)
    
    # PLOT ---
    p <- getBoxPlot(dm=dat)

    return(p)
}



###########################################################################################################################
### LOAD GENE ANNOTATION DATABASE FROM GFF/GTF ---
genecode.txdb <- GenomicFeatures::makeTxDbFromGFF(file=file.annotation_gtf, format="gtf", organism="Homo sapiens", chrominfo=seqinfo(TxDb.Hsapiens.UCSC.hg38.knownGene))


###########################################################################################################################
### LOAD TF HITS DATA ---
list.tf_hits <- readRDS(file=file.tf_hits_rds)
dat_tfsubtypes <- list.tf_hits$tf_footprints_hits_refined

### GET TF BY SUBTYPE ---
list.motif_subtype <- list()
for(subtype in subtypes){
    list.motif_subtype[[subtype]] <- get_motif_by_subtype(dat_tfsubtypes, subtype)
}   

### GET TF-PEAKS DIST TSS PLOT ----
list.plot <- list()
for(subtype in subtypes){
    tfs <- list.motif_subtype[[subtype]]
    list.plot[[subtype]] <- get_tf_peaks_distTSS(genecode.txdb, dir.chipatlas, chromosomes, tfs)
    cat("PLOT GENERATED:", subtype, "\n", sep=" ")
}





### GENERATE PLOT ---
file.plot <- file.path(dir.reproduce_fig, "supplementary_figure_18_c.pdf")
pdf(file.plot, width=2.5, height=4)
    grid.arrange(list.plot[[ subtypes[1] ]], nrow=1, ncol=1)
dev.off()

### GENERATE PLOT ---
file.plot <- file.path(dir.reproduce_fig, "supplementary_figure_18_d.pdf")
pdf(file.plot, width=2.5, height=1)
    grid.arrange(list.plot[[ subtypes[2] ]], nrow=1, ncol=1)
dev.off()

### GENERATE PLOT ---
file.plot <- file.path(dir.reproduce_fig, "supplementary_figure_18_e.pdf")
pdf(file.plot, width=2.5, height=2)
    grid.arrange(list.plot[[ subtypes[3] ]], nrow=1, ncol=1)
dev.off()

### GENERATE PLOT ---
file.plot <- file.path(dir.reproduce_fig, "supplementary_figure_18_f.pdf")
pdf(file.plot, width=2.5, height=4)
    grid.arrange(list.plot[[ subtypes[4] ]], nrow=1, ncol=1)
dev.off()

