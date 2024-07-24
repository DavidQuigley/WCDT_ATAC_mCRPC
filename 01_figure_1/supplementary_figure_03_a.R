###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### DEFINE LIBRARIES ---
library("stringr")
library("data.table")
library("dplyr")
library("ggplot2")
library("gridExtra")

### DEFINE PATH ---
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

### DEFINE FILES ----
file.rds_atacseq_masterdata <- file.path(dir.reproduce_data, "wcdt_atacseq_masterdata.rds")
file.utility_functions <- file.path(dir.reproduce_scripts, "00_prepare_data/00_01_utility_functions.R")
file.annotation_gtf <- file.path(dir.wrk, "reference/annotation_data/gene_annotation/gencode_v28/gencode.v28.annotation.gtf.gz")

### LOAD FUNCTIONS ---
source(file.utility_functions)

### LOAD ATAC-SEQ MASTER DATA ---
list.rds_atacseq_masterdata <- readRDS(file=file.rds_atacseq_masterdata)

### GET METADATA ---
metadata <- subset(list.rds_atacseq_masterdata$metadata_mcrpc, select=c("Sample_ID","ATACseq_Mapped_Reads","ATACseq_Peak_Counts","ATACseq_TSS_Enrich_Score","ATACseq_FRiP"))
colnames(metadata) <- c("SampleID","MappedReads","PeakCounts","TSSEnrichScore","FRiP")
metadata <- metadata[order(metadata$PeakCounts, decreasing=TRUE),]
rownames(metadata) <- NULL

### PREPARE DATA ---
dm <- reshape2::melt(metadata, id.vars="SampleID")
dm$SampleID <- factor(dm$SampleID, levels=metadata$SampleID)
dm$variable <- factor(dm$variable, levels=c("MappedReads","PeakCounts","TSSEnrichScore","FRiP"))


### GET PLOT ---
p1 <- ggplot(dm, aes(x=SampleID, y=value)) + 
            geom_bar(stat="identity", color=NA, fill="#006094",  width=0.8, size=0.5) + 
            facet_wrap(~variable, nrow=4, scales="free_y", drop=TRUE) +
            theme(
                axis.text.x = element_text(size = 6, color="#000000", angle=90, hjust=1, vjust=0.5),
                axis.text.y = element_text(size = 6, color="#000000"),
                axis.title = element_text(size = 10, color="#000000"),
                plot.title = element_text(size = 10, color="#000000", hjust=0.5),
                panel.grid.major.y = element_line(size=0.1, color="#BDBDBD"),
                panel.grid.major.x = element_line(size=0.1, color="#BDBDBD"),
                panel.grid.minor.x = element_blank(),
                panel.grid.minor.y = element_blank(),
                axis.ticks = element_line(size=0.2, color="#000000"), 
                strip.text = element_text(size=10, color="#000000"),
                strip.background = element_rect(fill="#FFFFFF", color="#FFFFFF"),
                panel.background = element_rect(fill="#FFFFFF", color="#000000"),
                legend.text = element_text(size = 10, color="#000000"),
                legend.title = element_blank(),
                legend.key.size = unit(0.5, "cm"),
                legend.position = "none") +
            ylab("Count") +
            xlab("") + 
            ggtitle("") 




######### PEAK ANNOTATION ##########################################################################################################

### FUNCTION: annotPeaksbyRange() ---
annotPeaksbyRange <- function(gr_peaks, genecode.txdb, sampleid){
    # ANNOTATE PEAKS ---
    annot_peaks <- ChIPseeker::annotatePeak(peak=gr_peaks, TxDb=genecode.txdb, verbose=TRUE)
    annot_peaks <- as.data.frame(annot_peaks)

    # ADD FEATURES ---
    annot_peaks$Feature <- ""
    annot_peaks$Feature[which(annot_peaks$annotation == "Promoter (<=1kb)")] <- "Promoter"
    annot_peaks$Feature[which(annot_peaks$annotation == "Promoter (1-2kb)")] <- "Promoter"
    annot_peaks$Feature[which(annot_peaks$annotation == "Promoter (2-3kb)")] <- "Promoter"
    annot_peaks$Feature[which(annot_peaks$annotation == "5' UTR")] <- "5' UTR"
    annot_peaks$Feature[which(annot_peaks$annotation == "3' UTR")] <- "3' UTR"

    exon <- annot_peaks$annotation[which(stringr::str_detect(annot_peaks$annotation, "Exon") == TRUE)]
    annot_peaks$Feature[which(annot_peaks$annotation %in% exon[which(str_detect(exon, "exon 1 of") == TRUE)])] <- "Exon"
    annot_peaks$Feature[which(annot_peaks$annotation %in% exon[which(str_detect(exon, "exon 1 of") == FALSE)])] <- "Exon"

    intn <- annot_peaks$annotation[which(stringr::str_detect(annot_peaks$annotation, "Intron") == TRUE)]
    annot_peaks$Feature[which(annot_peaks$annotation %in% intn[which(str_detect(intn, "intron 1 of") == TRUE)])] <- "Intron"
    annot_peaks$Feature[which(annot_peaks$annotation %in% intn[which(str_detect(intn, "intron 1 of") == FALSE)])] <- "Intron"

    annot_peaks$Feature[which(stringr::str_detect(annot_peaks$annotation, "Downstream") == TRUE)] <- "Downstream (<=300)"
    annot_peaks$Feature[which(annot_peaks$annotation == "Distal Intergenic")] <- "Distal Intergenic"

    # ADD SAMPLEID --
    annot_peaks$SampleID <- sampleid

    return(annot_peaks)
}


### LOAD GENE ANNOTATION DATABASE FROM GFF/GTF ---
genecode.txdb <- GenomicFeatures::makeTxDbFromGFF(file=file.annotation_gtf, format="gtf", organism="Homo sapiens", chrominfo=seqinfo(TxDb.Hsapiens.UCSC.hg38.knownGene))

### LOOP FOR EACH SAMPLEID ---
list.annot_peaks <- list()
for(i in 1:nrow(list.rds_atacseq_masterdata$metadata_mcrpc)){
    sample_name <- list.rds_atacseq_masterdata$metadata_mcrpc$Sample_Name[i]

    cat(format(Sys.time(), "%a %b %d %X %Y"), "ANNOTATION START:", sample_name, "\n", sep=" ")

    # GET PEAK BED ---
    dat_peaks <- list.rds_atacseq_masterdata$atacseq_peaks[[sample_name]]
    dat_peaks <- subset(dat_peaks, select=c("chr","start","end"))
    dat_peaks <- dat_peaks[!duplicated(dat_peaks),]

    # CONVERT PEAK BED TO GRANGES OBJECT ---
    gr_peaks <- GenomicRanges::makeGRangesFromDataFrame(df=dat_peaks, keep.extra.columns=FALSE)

    # GET PEAK ANNOTATION ---
    sampleid <- list.rds_atacseq_masterdata$metadata_mcrpc$Sample_ID[i]
    list.annot_peaks[[sampleid]] <- annotPeaksbyRange(gr_peaks, genecode.txdb, sampleid=sampleid)

    cat(format(Sys.time(), "%a %b %d %X %Y"), "DONE:", sample_name, "\n", sep=" ")
}

### AGGREGATE DATA ----
df.annot_peaks <- do.call(rbind.data.frame, list.annot_peaks)
rownames(df.annot_peaks) <- NULL

### GET FEATURE COUNT ---
dm <- df.annot_peaks %>%
        dplyr::select(SampleID, Feature) %>%
        dplyr::group_by(SampleID) %>%
        dplyr::count(Feature, sort=FALSE, name="Freq") %>%
        dplyr::mutate(Percentage=Freq/sum(Freq)*100)



### FEATURE ITEMS ---
feature_list <- c("Promoter","Distal Intergenic",
                    "Intron","Exon","5' UTR","3' UTR",
                    "Downstream (<=300)")

### COLOR ---
cpalette <- c("#e31a1c","#b15928",
                "#33a02c","#a6cee3","#ff7f00","#fdbf6f",
                "#cab2d6")

### FACTORIZE DATA ---
dm$SampleID <- factor(dm$SampleID, levels=metadata$SampleID)
dm$Feature <- factor(dm$Feature, levels=rev(feature_list))


### PLOT ---
p2 <- ggplot(dm, aes(x=SampleID, y=Percentage, fill=Feature)) + 
		geom_bar(stat="identity", color="#000000", width=0.8, size=0.2) + 
		scale_fill_manual(values=rev(cpalette)) +
		theme(
			axis.text.x = element_text(size = 6, color="#000000", angle=90, hjust=1, vjust=0.5),
			axis.text.y = element_text(size = 6, color="#000000"),
			axis.title = element_text(size = 8, color="#000000"),
			plot.title = element_text(size = 10, color="#000000", hjust=0),
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			axis.ticks = element_line(size = 0.2, colour="#000000"),	
			panel.background = element_rect(fill = "#FFFFFF", colour = "#000000"),           
			legend.text = element_text(size = 5, color="#000000"),
			legend.title = element_blank(),
			legend.key.size = unit(0.3, "cm"),			
			legend.position="bottom") + 
			ylab("Percentage") +
			xlab("") + 
            ggtitle("") 



################################################################################
### GENERATE PLOT ---
file.plot <- file.path(dir.reproduce_fig, "supplementary_figure_03_a.pdf")
pdf(file.plot, width=6.5, height=8)
    grid.arrange(p1, p2, nrow=2, ncol=1, heights=c(0.7,0.3))
dev.off()
