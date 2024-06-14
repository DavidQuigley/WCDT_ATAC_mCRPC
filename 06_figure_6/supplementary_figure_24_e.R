###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### LOAD LIBRARIES ---
library("stringr")
library("data.table")
library("GenomicRanges")
library("ggplot2")
library("gridExtra")

### DEFINE PATH ---
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/analysis/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

dir.footprints <- file.path(dir.reproduce_data, "tf_footprints/mcrpc_subtypes")

### DEFINE FILES ---
file.rds_footprints <- file.path(dir.footprints, "data_footprints/tobias_mcrpc_subtypes_tf_footprints_featureids_bound_unbound.rds")
file.rds_chipseq_bed <- file.path(dir.reproduce_data, "chipseq_bed.rds")

### MOTIF ---
motif_prefix_znf <- "MA0528.1_ZNF263"
motif_prefix_myc <- "MA0147.3_MYC"
motif_prefix_ar <- "ExtendedSite_AR"

motif_prefix <- motif_prefix_znf

### FUNCTION: getFootprintRanges() ---
getFootprintRanges <- function(feature_ids){
    # GET BED FILE ---
    df.bed <- data.frame(chr=unlist(lapply(stringr::str_split(feature_ids, "_"), function(x) x[1])),
                            start=as.numeric(unlist(lapply(stringr::str_split(feature_ids, "_"), function(x) x[2]))),
                            end=as.numeric(unlist(lapply(stringr::str_split(feature_ids, "_"), function(x) x[3]))))

    # REMOVE DUPLICATES ---
    df.bed <- df.bed[!duplicated(df.bed),]

    # GET GRANGES ---
    gr <- GenomicRanges::makeGRangesFromDataFrame(df.bed, keep.extra.columns=FALSE)

    return(gr)
}

### FUNCTION: callFootprint2ChipseqOverlap() ---
callFootprint2ChipseqOverlap <- function(list.chipseq, gr_footprint, analysis){
    list.dat <- list()
    for(i in 1:length(list.chipseq)){
        id <- names(list.chipseq)[i]
        p_overlap <- getFootprintOverlap(gr_footprint=gr_footprint, gr_chipseq=list.chipseq[[id]])
        list.dat[[id]] <- data.frame(dataset=id, analysis=analysis, value=p_overlap)
        cat("PROCESSED:", id, "\n", sep=" ")
    }

    # AGGREGATE DATA ---
    dat <- do.call(rbind.data.frame, list.dat)
    rownames(dat) <- NULL

    return(dat)
}

### FUNCTION: getFpChipSeqOverlap() ---
getFootprintOverlap <- function(gr_chipseq, gr_footprint){
    # GET OVERLAPS ---
    gr_overlap <- IRanges::subsetByOverlaps(x=gr_footprint, ranges=gr_chipseq)
    p_overlap <- round((length(gr_overlap) / length(gr_footprint) * 100), 2)

    return(p_overlap)
}

### LOAD CHIPSEQ DATA ---
list.rds_chipseq_bed <- readRDS(file=file.rds_chipseq_bed)
list.chipseq <- list.rds_chipseq_bed$chipseq_myc

### LOAD FOOTPRINTS ---
list.rds_footprints <- readRDS(file=file.rds_footprints)



### GET FOOTPRINT BED ---
gr_footprint_znf_bound <- getFootprintRanges(feature_ids=unique( unlist(list.rds_footprints$bound[[motif_prefix]]) ) )
gr_footprint_znf_unbound <- getFootprintRanges(feature_ids=unique( unlist(list.rds_footprints$unbound[[motif_prefix]]) ) )


### GET ALL BINDING SITES ---
gr_all <- c(gr_footprint_znf_bound, gr_footprint_znf_unbound)
gr_all <- gr_all[!(duplicated(gr_all))]



### GET OVERLAP OF ZNF263 FOOTPRINT TO AR CHIPSEQ  ---
df_atacseq_all <- callFootprint2ChipseqOverlap(list.chipseq=list.chipseq, gr_footprint=gr_all, analysis="all")
df_atacseq_bound <- callFootprint2ChipseqOverlap(list.chipseq=list.chipseq, gr_footprint=gr_footprint_znf_bound, analysis="bound")
df_atacseq_unbound <- callFootprint2ChipseqOverlap(list.chipseq=list.chipseq, gr_footprint=gr_footprint_znf_unbound, analysis="unbound")

### MERGE DATA ---
dm_atacseq <- rbind(df_atacseq_all, df_atacseq_bound, df_atacseq_unbound)

### FACTORIZE --
dm_atacseq$dataset <- factor(dm_atacseq$dataset, levels=rev( df_atacseq_bound$dataset[order(df_atacseq_bound$value, decreasing=TRUE)] ))
dm_atacseq$analysis <- factor(dm_atacseq$analysis, levels=rev(c("bound","unbound","all")))


### PLOT BAR ----
p <- ggplot(dm_atacseq, aes(x=value, y=dataset)) +
            geom_bar(aes(fill=analysis), stat="identity", position="dodge", color="#000000", width=0.8, size=0.2) +
            coord_cartesian(xlim=c(0,20)) +
            scale_x_continuous(breaks=seq(0,20, by=5) ) +
            scale_fill_manual(values=c("#969696","#a6cee3","#e31a1c")) +
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
               legend.key.size = unit(0.3, "cm"),
               legend.position = "bottom") +
            xlab("% of overlap of ZNF263 footprints (ATAC-seq)\n with that of MYC ChIP-seq peaks") +
            ylab("") + 
            ggtitle("")


### GENERATE PLOT ---
file.plot <- file.path(dir.reproduce_fig, "supplementary_figure_24_e.pdf")
pdf(file.plot, height=2, width=2)
    grid.arrange(p, nrow=1, ncol=1)
dev.off()

