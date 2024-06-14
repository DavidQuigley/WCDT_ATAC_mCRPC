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

dir.chipatlas <- file.path(dir.wrk, "external_data/chipatlas/data") # PATH CONTAINING CHIPATLAS FILES

### DEFINE FILES ---
file.annotation_gtf <- file.path(dir.wrk, "reference/annotation_data/gene_annotation/gencode_v28/gencode.v28.annotation.gtf.gz")

### CHROMOSOMES ---
chromosomes <- paste("chr", c(1:22,"X","Y"), sep="")


######################################################################################################
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


######################################################################################################
### TRANSCRIPTION FACTORS ---
tfs <- c("SP1","SP2","ZNF263","KLF5","SP3","MZF1","KLF16","CTCFL","KLF9","E2F6","SP4","RREB1","KLF14","PLAG1","NRF1","KLF4","ZNF740","EGR1","EGR2","GLIS2","NFKB2")

### DEFINE FILES ---
files.chipseq <- sapply(tfs, function(x) { paste(dir.chipatlas, "/Oth.ALL.05.", sprintf("%s", x), ".AllCell.bed", sep="") })

### GET CHIPSEQ PEAKS FOR EACH TF ---
list.chipseq <- list()
for(tf in tfs){
    cat("START:", tf, "\n", sep=" ")

    list.chipseq[[tf]] <- parse_chipseq_data(file.chipseq=files.chipseq[tf], chromosomes=chromosomes)

    cat("PROCESSED:", tf, "\n", sep=" ")
}

### CONVERT TO GRANGE LIST ---
list.gr_chipseq <- lapply(list.chipseq, function(x) GenomicRanges::makeGRangesFromDataFrame(x, keep.extra.columns=FALSE) )



### LOAD GENE ANNOTATION DATABASE FROM GFF/GTF ---
genecode.txdb <- GenomicFeatures::makeTxDbFromGFF(file=file.annotation_gtf, format="gtf", organism="Homo sapiens", chrominfo=seqinfo(TxDb.Hsapiens.UCSC.hg38.knownGene))

### GET MOTIF POSITION ANNOTATION ---
list.annot <- list()
for(motif_prefix in tfs){
    list.annot[[motif_prefix]] <- annotPeaksbyRange(gr_peaks=list.gr_chipseq[[motif_prefix]], genecode.txdb, motif_prefix)
    cat("PROCESSED", motif_prefix, "\n", sep=" ")
}


###########################################################################################################################
### AGGREGATE DATA ---
dat <- do.call(rbind.data.frame, list.annot)
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
        dplyr::arrange(Median)


### FACTORIZE ---
dat$motif_prefix <- factor(dat$motif_prefix, levels=df$motif_prefix)

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

# PLOT ---
p1 <- getBoxPlot(dm=dat)


###########################################################################################################################


### GET FEATURE COUNT ---
dm <- dat %>%
        dplyr::select(motif_prefix, Feature) %>%
        dplyr::group_by(motif_prefix) %>%
        dplyr::count(Feature, sort=FALSE, name="Freq") %>%
        dplyr::mutate(Percentage=Freq/sum(Freq)*100)

### FEATURE ITEMS ---
feature_list <- c("Promoter","Intergenic", "Intron","Exon","5UTR","3UTR","Downstream")

### COLOR ---
cpalette <- c("#e31a1c","#b15928","#33a02c","#a6cee3","#ff7f00","#fdbf6f","#cab2d6")

### FACTORIZE DATA ---
dm$motif_prefix <- factor(dm$motif_prefix, levels=df$motif_prefix)
dm$Feature <- factor(dm$Feature, levels=rev(feature_list))



### PLOT ---
p2 <- ggplot(dm, aes(x=Percentage, y=motif_prefix, fill=Feature)) + 
		geom_bar(stat="identity", color="#000000", width=0.8, size=0.2) + 
		scale_fill_manual(values=rev(cpalette)) +
		theme(
			axis.text.x = element_text(size = 6, color="#000000"),
			axis.text.y = element_text(size = 6, color="#000000"),
			axis.title = element_text(size = 8, color="#000000"),
			plot.title = element_text(size = 8, color="#000000", hjust=0.5),
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			axis.ticks = element_line(size = 0.2, colour="#000000"),	
			panel.background = element_rect(fill = "#FFFFFF", colour = "#000000"),           
			legend.text = element_text(size = 5, color="#000000"),
			legend.title = element_blank(),
			legend.key.size = unit(0.3, "cm"),			
			legend.position="none") + 
			ylab("") +
			xlab("Percentage") + 
            ggtitle("") 




### GENERATE PLOT ---
file.plot <- file.path(dir.reproduce_fig, "supplementary_figure_18_a_b.pdf")
pdf(file.plot, width=6, height=3)
    grid.arrange(p1, p2, nrow=1, ncol=2)
dev.off()

