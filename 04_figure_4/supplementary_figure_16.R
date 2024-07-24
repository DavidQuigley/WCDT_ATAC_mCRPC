###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### LOAD LIBRARIES ---
library("stringr")
library("data.table")
library("GenomicRanges")
library("dplyr")
library("ggplot2")
library("gridExtra")

### DEFINE PATH ---
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

dir.footprints <- file.path(dir.reproduce_data, "tf_footprints/mcrpc_subtypes")
dir.chipatlas <- file.path(dir.wrk, "external_data/chipatlas/data")

### DEFINE FILES ---
file.rds_footprints <- file.path(dir.footprints, "data_footprints/tobias_mcrpc_subtypes_tf_footprints_featureids_bound_unbound.rds")
file.rds_tf_hits <- file.path(dir.footprints, "data_footprints/tobias_mcrpc_subtypes_tf_hits.rds")


############################################################################################################
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

### FUNCTION: parse_chipseq_data() ---
parse_chipseq_data <- function(file.chipseq){
    # LOAD DATA ---
    dat <- data.table::fread(file=file.chipseq, sep="\t", header=FALSE, nThread=50, data.table=FALSE, verbose=FALSE)
    
    # ADD COLUMN NAMES ---
    colnames(dat)[1] <- "chr"
    colnames(dat)[2] <- "start"
    colnames(dat)[3] <- "end"

    return(dat)
}

### FUNCTION: getFootprintOverlap() ---
getFootprintOverlap <- function(gr_chipseq, gr_footprint){
    # GET OVERLAPS ---
    gr_overlap <- IRanges::subsetByOverlaps(x=gr_footprint, ranges=gr_chipseq)
    p_overlap <- (length(gr_overlap) / length(gr_footprint) * 100)

    return(p_overlap)
}

### FUNCTION: getPlot() ---
getPlot <- function(dm){
    # PLOT ---
    p <- ggplot(dm, aes(x=value, y=gene)) +
            #geom_bar(aes(fill=analysis), stat="identity", position="dodge", color="#000000", width=0.8, size=0.2) +
            geom_bar(fill="#e31a1c", stat="identity", position="dodge", color="#000000", width=0.8, size=0.2) +
            coord_cartesian(xlim=c(0,100)) +
            scale_x_continuous(breaks=seq(0,100, by=10) ) +
            #scale_fill_manual(values=c("#a6cee3","#e31a1c")) +
            theme(
               axis.text.x = element_text(size = 6, color="#000000"),
               axis.text.y = element_text(size = 6, color="#000000"),
               axis.title = element_text(size = 8, color="#000000"),
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
               legend.position = "none") +
            xlab("% of overlap of TF footprints (ATAC-seq)\n with that of ChIP-seq peaks") +
            ylab("") + 
            ggtitle("")

    return(p)
}

############################################################################################################
fnames <- list.files(dir.chipatlas, pattern=".bed", full.names=FALSE)
tfs <- unlist(lapply(stringr::str_split(fnames, "[.]"), function(x) x[4]))

############################################################################################################
### LOAD TF MOTIF DATA ---
list.tf_hits <- readRDS(file=file.rds_tf_hits)
des <- list.tf_hits$tf_footprints_hits_refined
des <- subset(des, select=c("motif_prefix","gene"))
des$gene[which(des$gene == "Pax6")] <- "PAX6"
des <- subset(des, des$gene %in% tfs)

### DEFINE CHIPSEQ FILES ---
files.chipseq <- sapply(tfs, function(x) { paste(dir.chipatlas, "/Oth.ALL.05.", sprintf("%s", x), ".AllCell.bed", sep="") })

### LOAD FOOTPRINTS ---
list.rds_footprints <- readRDS(file=file.rds_footprints)



############################################################################################################
### GET BOUND FOOTPRINTS GRANGES ---
list.gr_motif_bound <- list()
for(motif_prefix in des$motif_prefix){
    # GET MOTIF FEATURES ---
    feature_ids <- unique( unlist(list.rds_footprints$bound[[motif_prefix]]) )

    # GET MOTIF RANGES ---
    list.gr_motif_bound[[motif_prefix]] <- getGRanges(feature_ids)
}

### GET UNBOUND FOOTPRINTS GRANGES ---
list.gr_motif_unbound <- list()
for(motif_prefix in des$motif_prefix){
    # GET MOTIF FEATURES ---
    feature_ids <- unique( unlist(list.rds_footprints$unbound[[motif_prefix]]) )

    # GET MOTIF RANGES ---
    list.gr_motif_unbound[[motif_prefix]] <- getGRanges(feature_ids)
}



############################################################################################################
### LOOP FOR EACH TF ---
list.dm <- list()
for(i in 1:nrow(des)){
    motif_prefix <- des$motif_prefix[i]
    gene <- des$gene[i]
    
    cat("START:", gene, "\n", sep=" ")

    # LOAD CHIP-SEQ DATA ---
    dat_chipseq <- parse_chipseq_data(file.chipseq=files.chipseq[gene])
    gr_chipseq <- GenomicRanges::makeGRangesFromDataFrame(dat_chipseq, keep.extra.columns=FALSE)

    # ATAC-SEQ TF FOOTPRINTS ---
    gr_motif_bound <- list.gr_motif_bound[[motif_prefix]]
    gr_motif_unbound <- list.gr_motif_unbound[[motif_prefix]]

    p_overlap_bound <- getFootprintOverlap(gr_chipseq=gr_chipseq, gr_footprint=gr_motif_bound)
    p_overlap_unbound <- getFootprintOverlap(gr_chipseq=gr_chipseq, gr_footprint=gr_motif_unbound)

    list.dm[[motif_prefix]] <- data.frame(motif_prefix=motif_prefix, gene=gene, analysis=c("bound","unbound"), value=c(p_overlap_bound,p_overlap_unbound))
    
    cat("END:", gene, "\n", sep=" ")
}

### AGGREGATE DATA ---
dm <- do.call(rbind.data.frame, list.dm)
rownames(dm) <- NULL

### WRITE OUTPUT ---
#file.output <- file.path(dir.wrk, "external_data/chipatlas/parsed/chipatlas_atacseq_footprint_comparison.tsv")
#write.table(dm, file.output, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

### RE-LOAD DATA --
#file.dm <- file.path(dir.wrk, "external_data/chipatlas/parsed/chipatlas_atacseq_footprint_comparison.tsv")
#dm <- data.table::fread(file=file.dm, sep="\t", header=TRUE, nThread=1, data.table=FALSE, verbose=FALSE)


### RETAIN RESULTS FOR TFs WITH MORE THAN 5% OF ATAC-SEQ TF FOOTPRINT OVERLAP WITH CHIP-SEQ DATA ----
tf_pass <- dm$gene[which( (dm$analysis == "bound") & (dm$value > 5))]
dm_pass <- subset(dm, dm$gene %in% tf_pass)
dm_pass <- subset(dm_pass, dm_pass$analysis == "bound")

### ORDER GENES ---
d <- subset(dm_pass, dm_pass$analysis == "bound")
d <- d[order(d$value, decreasing=TRUE),]

d1 <- d[1:60,]
d2 <- d[61:120,]

### PARTITION DATA ---
dm_pass_1 <- subset(dm_pass, dm_pass$gene %in% d1$gene)
dm_pass_2 <- subset(dm_pass, dm_pass$gene %in% d2$gene)

### FACTORIZE ---
dm_pass_1$gene <- factor(dm_pass_1$gene, levels=rev(d1$gene))
#dm_pass_1$analysis <- factor(dm_pass_1$analysis, levels=rev(c("bound","unbound")))

dm_pass_2$gene <- factor(dm_pass_2$gene, levels=rev(d2$gene))
#dm_pass_2$analysis <- factor(dm_pass_2$analysis, levels=rev(c("bound","unbound")))


# PLOT ---
p1 <- getPlot(dm=dm_pass_1)
p2 <- getPlot(dm=dm_pass_2)

### GENERATE PLOT ---
file.plot <- file.path(dir.reproduce_fig, "supplementary_figure_16.pdf")
pdf(file.plot, height=6, width=6)
    grid.arrange(p1, p2, nrow=1, ncol=2)
dev.off()
