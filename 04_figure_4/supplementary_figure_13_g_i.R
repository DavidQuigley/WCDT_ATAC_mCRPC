### LOAD LIBRARIES ---
library("stringr")
library("data.table")
library("GenomicRanges")
library("IRanges")
library("ggplot2")
library("gridExtra")

### DEFINE PATH ---
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

dir.footprints <- file.path(dir.reproduce_data, "tf_footprints/mcrpc_subtypes")
dir.bam <- file.path(dir.reproduce_data, "tf_footprints/mcrpc_subtypes/bam")
dir.bed <- file.path(dir.reproduce_data, "tf_footprints/mcrpc_subtypes/bed")

### DEFINE FILES ---
file.rds_footprints <- file.path(dir.footprints, "data_footprints/expt_bound_pvalue_footprint_summary.rds")
file.rds_chipseq <- file.path(dir.reproduce_data, "chipseq_bed.rds")

### SET VALUES ---
val_analysis <- "ARpNEn_ARnNEp"
val_pvalue <- "p001"
motifs <- c("ExtendedSite_AR","MA0148.3_FOXA1","MA0901.1_HOXB13")


### FUNCTION: loadFootprintData() ---
loadFootprintData <- function(list.rds_footprints, val_pvalue, val_analysis, motif){
    # LOAD DATA ---
    dat <- list.rds_footprints[[val_pvalue]][[val_analysis]][[motif]]
    
    # FIX COLUMN NAMES ---
    colnames(dat)[1] <- "chr"
    colnames(dat)[2] <- "start"
    colnames(dat)[3] <- "end"

    # GET BOUND/UNBOUND COLUMN NAME ---
    val_cname <- paste(stringr::str_split(val_analysis, "_")[[1]][1], "bound", sep="_")

    # SEPARATE BOUND AND UNBOUND DATA ---
    dat_bound <- dat[which(dat[,val_cname] == 1),]
    dat_unbound <- dat[which(dat[,val_cname] == 0),]

    # GET GRANGES ---
    gr_bound <- GenomicRanges::makeGRangesFromDataFrame(dat_bound, keep.extra.columns=TRUE)
    gr_unbound <- GenomicRanges::makeGRangesFromDataFrame(dat_unbound, keep.extra.columns=TRUE)

    # COMPILE INTO LIST ---
    list.gr <- list()
    list.gr$bound <- gr_bound
    list.gr$unbound <- gr_unbound

    return(list.gr)
}

### FUNCTION: get_overlap() ---
get_overlap <- function(gr.atacseq, gr.chipseq){
    # GET OVERLAP ---
    gr.overlap_atacseq <- IRanges::subsetByOverlaps(x=gr.atacseq, ranges=gr.chipseq, type="any")

    return(gr.overlap_atacseq)
}    



### LOAD ATAC-SEQ FOOTPRINTS DATA ---
list.rds_footprints <- readRDS(file=file.rds_footprints)

### LOAD CHIP-SEQ DATA ---
list.rds_chipseq <- readRDS(file=file.rds_chipseq)

### GET CHIP-SEQ DATA TO LIST ---
gr_chipseq_AR <- list.rds_chipseq$pomerantz_mcrpc$AR
gr_chipseq_FOXA1 <- list.rds_chipseq$pomerantz_mcrpc$FOXA1
gr_chipseq_HOXB13 <- list.rds_chipseq$pomerantz_mcrpc$HOXB13
list.chipseq <- list(gr_chipseq_AR, gr_chipseq_FOXA1, gr_chipseq_HOXB13)
names(list.chipseq) <- motifs




### GET TF-FOOTPRINT AND CHIP-SEQ OVERLAP ----------
list.res <- list()
for(motif in motifs){
    # GET CHIP-SEQ DATA ---
    gr_chipseq <- list.chipseq[[motif]]

    # GET BOUND/UNBOUND FOOTPRINT DATA ---
    list.gr_footprint <- loadFootprintData(list.rds_footprints, val_pvalue, val_analysis, motif)
    gr_bound <- list.gr_footprint$bound
    gr_unbound <- list.gr_footprint$unbound

    # GET OVERLAP ---
    gr_overlap_bound <- get_overlap(gr.atacseq=gr_bound, gr.chipseq=gr_chipseq)
    gr_overlap_unbound <- get_overlap(gr.atacseq=gr_unbound, gr.chipseq=gr_chipseq)

    # CALC STATS ---
    n_bound <- length(gr_bound)
    n_unbound <- length(gr_unbound)
    n_overlap_bound <- length(gr_overlap_bound) 
    n_overlap_unbound <- length(gr_overlap_unbound) 

    p_overlap_bound <- round( (n_overlap_bound/n_bound) * 100, 2)
    p_overlap_unbound <- round( (n_overlap_unbound/n_unbound) * 100, 2)

    # PREPARE DATA ---
    list.res[[motif]] <- rbind( data.frame(motif=motif, status="bound", p_overlap=p_overlap_bound),
                                data.frame(motif=motif, status="unbound", p_overlap=p_overlap_unbound) )


    cat("PROCESSED:", motif, "\n", sep=" ")
}

### AGGREGATE DATA ---
df <- do.call(rbind.data.frame, list.res)
rownames(df) <- NULL

### FACTORIZE --
df$motif <- factor(df$motif, levels=motifs )
df$status <- factor(df$status, levels=rev(c("bound","unbound")))


### PLOT BAR ----
p <- ggplot(df, aes(x=p_overlap, y=status)) +
            geom_bar(aes(fill=status), stat="identity", position="dodge", color="#000000", width=0.8, size=0.2) +
            coord_cartesian(xlim=c(0,100)) +
            scale_x_continuous(breaks=seq(0,100, by=10) ) +
            scale_fill_manual(values=c("#a6cee3","#e31a1c")) +
            facet_wrap(~motif, nrow=1, ncol=3, scales="fixed") + 
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
            xlab("% of overlap of TF footprints (ATAC-seq)\n with that of ChIP-seq peaks") +
            ylab("") + 
            ggtitle("")

### GENERATE PLOT ---
file.plot <- file.path(dir.reproduce_fig, "supplementary_figure_11_g_i.pdf")
pdf(file.plot, height=2, width=6)
    grid.arrange(p, nrow=1, ncol=1)
dev.off()





#> df
#            motif  status p_overlap
#1 ExtendedSite_AR   bound     82.17
#2 ExtendedSite_AR unbound     33.67
#3  MA0148.3_FOXA1   bound     96.97
#4  MA0148.3_FOXA1 unbound     65.20
#5 MA0901.1_HOXB13   bound     95.12
#6 MA0901.1_HOXB13 unbound     60.70

