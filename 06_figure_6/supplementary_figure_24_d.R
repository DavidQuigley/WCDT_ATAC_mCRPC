###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### LOAD LIBRARIES ---
library("stringr")
library("data.table")
library("reshape2")
library("IRanges")
library("GenomicRanges")
library("ggplot2")
library("gridExtra")

### DEFINE PATH ---
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

dir.tfpeaks2gene <- file.path(dir.reproduce_data, "tfpeaks2gene")

### DEFINE FILES ---
file.rds_tf2peaks <- file.path(dir.tfpeaks2gene, "wcdt_tf2peaks_mcrpc_subtypes.rds")
file.utility_functions <- file.path(dir.reproduce_scripts, "00_prepare_data/00_01_utility_functions.R")

### LOAD FUNCTIONS ---
source(file.utility_functions)

### MOTIFS ---
motif_prefix_znf <- "MA0528.1_ZNF263"
motif_prefix_myc <- "MA0147.3_MYC"

### CHROMOSOMES ---
chromosome <- paste("chr", c(1:22,"X","Y"), sep="")

### SUBTYPES ---
subtypes <- c("ARpNEn","ARlNEn","ARpNEp","ARnNEp","ARnNEn")

### COLOR ---
cp_subtypes <- c("#e31a1c","#fb9a99","#ff7f00","#33a02c","#1f78b4")

### BINS ---
bin_windows <- seq(1000,50, by=-50)


############################################################################################################################################
### FUNCTION: getGranges() ---
getGranges <- function(df){
    # GET CO-ORDINATES ---
    df$chr <- unlist(lapply(stringr::str_split(df$motif_position, "_"), function(x) x[1]))
    df$start <- as.numeric(unlist(lapply(stringr::str_split(df$motif_position, "_"), function(x) x[2])))
    df$end <- as.numeric(unlist(lapply(stringr::str_split(df$motif_position, "_"), function(x) x[3])))

    # GET GRANGES OBJECT ---
    gr <- GenomicRanges::makeGRangesFromDataFrame(df, keep.extra.columns=TRUE)

    return(gr)
}

### FUNCTION: getOverlapsRanges() ---
getOverlapsRanges <- function(df_znf, df_myc){
    # GET GRANGES OBJECT ---
    gr_znf <- getGranges(df=df_znf)
    gr_myc <- getGranges(df=df_myc)

    # GET OVERLAPS ---
    gr_overlap <- IRanges::subsetByOverlaps(x=gr_znf, ranges=gr_myc)

    return(gr_overlap)
}

### FUNCTION: getFootprintRanges() ---
getFootprintRanges <- function(list.tf2peaks, motif_prefix, subtype){
    # GET DATA ---
    dat <- list.tf2peaks[[motif_prefix]][[subtype]]

    # GET GENOMEIC COORDINATES: TF MOTIF ---
    dat$chr <- unlist(lapply(stringr::str_split(dat$motif_position, "_"), function(x) x[1]))
    dat$start <- as.numeric(unlist(lapply(stringr::str_split(dat$motif_position, "_"), function(x) x[2])))
    dat$end <- as.numeric(unlist(lapply(stringr::str_split(dat$motif_position, "_"), function(x) x[3])))

    # GET GRANGES ---
    gr_motif <- GenomicRanges::makeGRangesFromDataFrame(df=dat, keep.extra.columns=FALSE)

    return(gr_motif)
}

### FUNCTION: getMotifBins() ---
getMotifBins <- function(bins, list.tf2peaks, motif_prefix, subtype){
    # GET GRANGES ---
    gr_motif <- getFootprintRanges(list.tf2peaks, motif_prefix, subtype)

    # GET OVERLAP OF BINS TO MOTIFS ---
    gr_overlap <- IRanges::mergeByOverlaps(query=bins, subject=gr_motif)
    bin_features <- paste(as.character(gr_overlap$bins@seqnames), as.character(gr_overlap$bins@ranges), sep="_")

    return(bin_features)
}

### FUNCTION: getTFbinOverlap() ---
getTFbinOverlap <- function(bins, bin_window, list.tf2peaks, list.overlap_ranges_znf_myc, motif_prefix_myc, motif_prefix_znf, subtype){
    # GET BIN FEATURES: MYC ---
    bin_features_myc <- getMotifBins(bins, list.tf2peaks, motif_prefix=motif_prefix_myc, subtype)

    # GET ZNF263-MYC OVERLAPING FOOTPRING REGION ---
    gr_znf_myc_overlap <- list.overlap_ranges_znf_myc[[subtype]]

    # GET BIN FEATURES: CO-FACTORS ---
    gr_znf <- getFootprintRanges(list.tf2peaks, motif_prefix=motif_prefix_znf, subtype)
    gr_znf_nonoverlap <- subsetByOverlaps(x=gr_znf, ranges=gr_znf_myc_overlap, invert=TRUE)

    gr_znf_bin <- IRanges::mergeByOverlaps(query=bins, subject=gr_znf_nonoverlap)
    bin_features_cofactors <- paste(as.character(gr_znf_bin$bins@seqnames), as.character(gr_znf_bin$bins@ranges), sep="_")


    # GET INTERSECTION WITH MYC MOTIF BINS ---
    percent_overlap <- round( (length(intersect(bin_features_cofactors, bin_features_myc))/length(bin_features_myc)) * 100, 2 )


    # PREPARE DATA ---
    df <- data.frame(bin=bin_window, subtype=subtype, percent_overlap=percent_overlap)

    return(df)
}

### FUNCTION: getOverlapsByBinWindow() ---
getOverlapsByBinWindow <- function(bins, bin_window, list.tf2peaks, list.overlap_ranges_znf_myc, motif_prefix_myc, motif_prefix_znf, subtypes){
    # LOOP BY SUBTYPES ---
    list.binoverlaps <- list()
    for(k in 1:length(subtypes)){
        subtype <- subtypes[k]
        list.binoverlaps[[subtype]] <- getTFbinOverlap(bins, bin_window, list.tf2peaks, list.overlap_ranges_znf_myc,motif_prefix_myc, motif_prefix_znf, subtype)
        cat("PROCESSED:", subtype, "\n", sep=" ")
    }

    # AGGREGATE DATA ---
    df <- do.call(rbind.data.frame, list.binoverlaps)
    rownames(df) <- NULL

    return(df)
}


### FUNCTION: getplot() ---
getplot <- function(df, color_palette){
    # PLOT ---
    p <- ggplot(df, aes(x=bin, y=percent_overlap)) +
            geom_line(aes(color=subtype), size=0.5, alpha=0.7) +
            geom_point(aes(color=subtype, fill=subtype), shape=21, stroke=0.25, size=1, alpha=0.7, na.rm=TRUE) +
            coord_cartesian(xlim=c(0,1000), ylim=c(0,60) ) +
            scale_x_continuous(breaks=seq(0, 1000, by=100) ) +
            scale_y_continuous(breaks=seq(0, 60, by=10) ) +
            scale_fill_manual(values=color_palette) +
            scale_color_manual(values=color_palette) +
            theme(
               axis.text.x = element_text(size = 6, color="#000000"),
               axis.text.y = element_text(size = 6, color="#000000"),
               axis.title = element_text(size = 6, color="#000000"),
               plot.title = element_text(size = 8, color="#000000", hjust=0.5),
               panel.grid.major.y = element_blank(),
               panel.grid.major.x = element_blank(),
               panel.grid.minor = element_blank(),
               axis.ticks.y = element_line(size=0.2, color="#000000"),
               axis.ticks.x = element_line(size=0.2, color="#000000"),
               strip.text = element_text(size=7, color="#000000"),
               strip.background = element_rect(fill="#FFFFFF", color="#FFFFFF"),
               panel.background = element_rect(fill="#FFFFFF", color="#000000"),
               legend.text = element_text(size = 5, color="#000000"),
               legend.title = element_blank(),
               legend.key.size = unit(0.4, "cm"),
               legend.position = "none") +
            xlab("Genome Bin Width (bp)") +
            ylab("% overlap of ZNF263-MYC footprints regions") + 
            ggtitle("ZNF263-MYC Non-overlapping")

    return(p)            
}



############################################################################################################################################
### LOAD TF2PEAKS ---
list.tf2peaks <- readRDS(file=file.rds_tf2peaks)



############################################################################################################################################
### GET ZNF263-MYC FOOTPRINT OVERLAP ---
list.overlap_ranges_znf_myc <- list()
for(subtype in subtypes){
    df_znf <- list.tf2peaks[[motif_prefix_znf]][[subtype]]
    df_myc <- list.tf2peaks[[motif_prefix_myc]][[subtype]]

    # GET OVERLAPS LIST ---
    list.overlap_ranges_znf_myc[[subtype]] <- getOverlapsRanges(df_znf, df_myc)
}



### GET OVERLAP OF ZNF263:MYC FOOTPRINTS BY BINS -----------
list.result <- list()
for(ctr in 1:length(bin_windows)){
    bin_window <- bin_windows[ctr]
    y <- paste("bin", bin_window, sep="_")

    cat("START BIN:", bin_window, "\n", sep=" ")

    # BIN GENOME ---
    bins <- binGenome(genome=BSgenome.Hsapiens.UCSC.hg38, chromosome=chromosome, tile.size=bin_window)

    # GET OVARLAPS BY BIN WINDOW ---
    list.result[[y]] <- getOverlapsByBinWindow(bins, bin_window, list.tf2peaks, list.overlap_ranges_znf_myc, motif_prefix_myc, motif_prefix_znf, subtypes)

    cat("PROCESSED BIN:", bin_window, "\n", sep=" ")
}

### AGGREGATE DATA ---
df <- do.call(rbind.data.frame, list.result)
rownames(df) <- NULL




############################################################################################################################################
### FACTORIZE ---
df$subtype <- factor(df$subtype, levels=subtypes)

### GET PLOT ---
p <- getplot(df, color_palette=cp_subtypes)

### GENERATE PLOT ---
file.plot <- file.path(dir.reproduce_fig, "supplementary_figure_24_d.pdf")
pdf(file.plot, height=2.5, width=2.5)
    grid.arrange(p, nrow=1, ncol=1)
dev.off()


