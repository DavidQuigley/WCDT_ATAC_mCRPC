###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### LOAD LIBRARIES ---
library("stringr")
library("data.table")
library("tibble")
library("dplyr")
library("ggplot2")
library("gridExtra")
library("ggrepel")

### DEFINE PATH ---
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/analysis/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

dir.tfpeaks2gene <- file.path(dir.reproduce_data, "tfpeaks2gene")
dir.footprints <- file.path(dir.reproduce_data, "tf_footprints/mcrpc_subtypes/data_footprints") 

### DEFINE FILES ---
file.tf2genes <- file.path(dir.tfpeaks2gene, "wcdt_tf2genes_r_0p4_pval_0p05_distTSS_500kb.tsv.gz")
file.tf_hits_rds <- file.path(dir.footprints, "tobias_mcrpc_subtypes_tf_hits.rds")

### SUBTYPES --
subtypes <- c("ARpNEn","ARlNEn","ARpNEp","ARnNEp","ARnNEn")




###########################################################################################################################
### FUNCTION: get_tf2gene_pairs() ---
get_tf2gene_pairs <- function(dat_tf2genes, dat_tfsubtypes, subtype){
    # EXTRACT TF BY SUBTYPE ---
    index_subtype <- which( colnames(dat_tfsubtypes) == subtype )
    index_tf <- which( dat_tfsubtypes[,index_subtype] >= 2 ) 
    dat_tfsubtypes <- dat_tfsubtypes[index_tf,]
    tfs_subtype <- dat_tfsubtypes$motif_prefix

    # EXTRACT TF2GENES BY TF ---
    dat_tf2genes_sub <- subset(dat_tf2genes, dat_tf2genes$motif_prefix %in% tfs_subtype)

    # EXTRACT TF-GENE PAIRS ---
    #df <- subset(dat_tf2genes_sub, select=c("motif_prefix","Gene"))
    #df <- df[!duplicated(df),]

    return(dat_tf2genes_sub)
}

### FUNCTION: prepareData() ---
prepareData <- function(df){
    # GET COUNT TABLE ---
    dm <- tibble::tibble(df) %>% 
            dplyr::select(motif_prefix, Gene) %>%
            dplyr::distinct(motif_prefix, Gene, .keep_all = TRUE) %>%
            dplyr::count(motif_prefix, sort=TRUE, name="Freq")

    dm$tf_name <- unlist(lapply(str_split(dm$motif_prefix, "_"), function(x) paste(x[2:length(x)], collapse="_") ))
    dm$Rank <- c(1:nrow(dm))    

    return(dm)
}


### FUNCTION: getRankPlot() ---
getRankPlot <- function(dm, list.params, subtype){
    dm$Gene_Label <- dm$tf_name

    # PARAMS ---
    x_low <- list.params[[subtype]]$x_low
    x_high <- list.params[[subtype]]$x_high
    x_interval <- list.params[[subtype]]$x_interval
    y_low <- list.params[[subtype]]$y_low
    y_high <- list.params[[subtype]]$y_high
    y_interval <- list.params[[subtype]]$y_interval

    # PLOT FREQ- RANK ----
    p1 <- ggplot(dm, aes(x=Rank, y=Freq, label=Gene_Label)) +
                geom_line(color="#006094", size=0.5, alpha=0.7) +
                geom_point(fill="#006094", color="#006094", shape=21, stroke=0.25, size=0.7, alpha=0.7, na.rm=TRUE) +
                coord_cartesian(xlim=c(x_low,x_high), ylim=c(y_low,y_high) ) +
                scale_x_continuous(breaks=seq(x_low, x_high, by=x_interval) ) +
                scale_y_continuous(breaks=seq(y_low, y_high, by=y_interval) ) +
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
                    strip.text = element_text(size=8, color="#000000"),
                    strip.background = element_rect(fill="#FFFFFF", color="#FFFFFF"),
                    panel.background = element_rect(fill="#FFFFFF", color="#000000"),
                    legend.text = element_text(size = 10, color="#000000"),
                    legend.title = element_blank(),
                    legend.key.size = unit(0.5, "cm"),
                    legend.position = "none") +
                xlab("Rank") +
                ylab("No. of target genes") + 
                ggtitle(subtype)

    p2 <- p1 + ggrepel::geom_label_repel(aes(label=Gene_Label), segment.color="#BDBDBD", color="#000000", 
                        alpha=0.8, size=1, label.size=NA, fill=NA, min.segment.length=0.001, max.overlaps=Inf, na.rm=TRUE)            

    # GATHER PLOTS TO LIST ---
    list.output <- list(p1=p1, p2=p2)

    return(list.output)
}



###########################################################################################################################
### LOAD TF2GENES DATA ---
dat_tf2genes <- data.table::fread(file=file.tf2genes, sep="\t", header=TRUE, nThread=50, data.table=FALSE, verbose=FALSE)

### LOAD TF HITS DATA ---
list.tf_hits <- readRDS(file=file.tf_hits_rds)
dat_tfsubtypes <- list.tf_hits$tf_footprints_hits_refined

### PARAMS ---
list.params <- list()
list.params[["ARpNEn"]] <- list(x_low=0, x_high=50, x_interval=10, y_low=0, y_high=3500, y_interval=500)
list.params[["ARlNEn"]] <- list(x_low=0, x_high=15, x_interval=5, y_low=0, y_high=5500, y_interval=500)
list.params[["ARpNEp"]] <- list(x_low=0, x_high=25, x_interval=5, y_low=0, y_high=3500, y_interval=500)
list.params[["ARnNEp"]] <- list(x_low=0, x_high=50, x_interval=10, y_low=0, y_high=4000, y_interval=500)
list.params[["ARnNEn"]] <- list(x_low=0, x_high=25, x_interval=5, y_low=0, y_high=5500, y_interval=500)




###########################################################################################################################
### GET TF2GENE PAIRS PER SUBTYPE ----
list.tf2gene <- list()
list.tf2gene_stats <- list()
list.plot <- list()
for(subtype in subtypes){
    list.tf2gene[[subtype]] <- get_tf2gene_pairs(dat_tf2genes, dat_tfsubtypes, subtype)
    list.tf2gene_stats[[subtype]] <- prepareData(df=list.tf2gene[[subtype]])
    list.plot[[subtype]] <- getRankPlot(dm=list.tf2gene_stats[[subtype]], list.params, subtype)

    cat("PROCESSED:", subtype, "\n", sep=" ")
}


###########################################################################################################################
### GENERATE PLOT ---
file.plot <- file.path(dir.reproduce_fig, "figure_5.pdf")
pdf(file.plot, height=1.7, width=1.7)
    for(subtype in subtypes){
        grid.arrange(list.plot[[subtype]]$p1, nrow=1, ncol=1)
        grid.arrange(list.plot[[subtype]]$p2, nrow=1, ncol=1)
    }
dev.off()
