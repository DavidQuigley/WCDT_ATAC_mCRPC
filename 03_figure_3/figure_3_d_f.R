###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### DEFINE LIBRARIES ---
library("stringr")
library("data.table")
library("GenomicRanges")
library("rGREAT")
library("ggplot2")
library("gridExtra")

### DEFINE PATH ---
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/analysis/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

### DEFINE FILES ----
file.rds_atacseq_read_counts_mcrpc <- file.path(dir.reproduce_data, "atacseq_read_counts_mcrpc.rds")
file.kwtest <- file.path(dir.reproduce_data, "wcdt_atacseq_mcrpc_subtypes_KruskalWallisTest_summary.tsv.gz")




### FUNCTION: getGRanges() ---
getGRanges <- function(dat){
    dat <- subset(dat, select=c("seqnames","start","end"))
    gr <- GenomicRanges::makeGRangesFromDataFrame(dat, keep.extra.columns=FALSE)
    return(gr)
}


### FUNCTION: processEnrichment() ---
processEnrichment <- function(dat, pval){
    dat <- subset(dat, select=c("id","fold_enrichment","p_value","p_adjust"))
    dat <- dat[which(dat$p_adjust <= pval),]
    dat <- dat[order(dat$p_value, decreasing=FALSE),]
    dat$nlogp <- -log10(dat$p_adjust)

    return(dat)
}

### FUNCTION: plotBar() ----
plotBar <- function(dat, analysis_id){
    # PREPARE DATA ---
    dat <- dat[order(dat$nlogp, decreasing=FALSE),]
    dat$id <- factor(dat$id, levels=dat$id)

    # PLOT ---
    p <- ggplot(dat, aes(x=nlogp, y=id)) +
            geom_bar(stat="identity", fill="#006094", width=0.8) +
            coord_cartesian(xlim=c(0,12)) +
            scale_x_continuous(breaks=seq(0,12,by=4)) +
            theme(
                axis.text.x = element_text(size = 5, color="#000000"),
                axis.text.y = element_text(size = 5, color="#000000"),
                axis.title = element_text(size = 5, color="#000000"),
                plot.title = element_text(size = 5, color="#000000", hjust=0.5),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.ticks = element_line(size=0.2, color="#000000"), 
                strip.text = element_text(size=5, color="#000000"),
                strip.background = element_rect(fill="#FFFFFF", color="#FFFFFF"),
                panel.background = element_rect(fill="#FFFFFF", color="#000000"),
                legend.text = element_text(size = 5, color="#000000"),
                legend.title = element_blank(),
                legend.key.size = unit(0.3, "cm"),
                legend.position = "none") +
            ylab("") +
            xlab("-log10(fdr)") + 
            ggtitle(analysis_id) 

    return(p)
}

#############################################################################################################################
### LOAD ATAC-SEQ PEAK ANNOTATION DATA ---
list.rds_atacseq_read_counts_mcrpc <- readRDS(file=file.rds_atacseq_read_counts_mcrpc)
annot <- list.rds_atacseq_read_counts_mcrpc$feature_annotation

### LOAD KRUSKAL-WALLIS TEST DATA ---
dat_kwtest <- data.table::fread(file=file.kwtest, sep="\t", header=TRUE, nThread=50, data.table=FALSE, verbose=FALSE)
colnames(dat_kwtest)[1] <- "FeatureID"
dat_kwtest <- merge(dat_kwtest, annot, by="FeatureID")
dat_kwtest <- dat_kwtest[order(dat_kwtest$pvalue, decreasing=FALSE),]
dat_kwtest <- subset(dat_kwtest, dat_kwtest$pvalue <= 0.001)
dat_kwtest <- subset(dat_kwtest, dat_kwtest$GeneType == "protein_coding")

dat_kwtest <- subset(dat_kwtest, dat_kwtest$Feature %in% c("Promoter","Intron","Distal Intergenic"))
dat_kwtest$Feature[which(dat_kwtest$Feature == "Promoter")] <- "promoter"
dat_kwtest$Feature[which(dat_kwtest$Feature == "Intron")] <- "intron"
dat_kwtest$Feature[which(dat_kwtest$Feature == "Distal Intergenic")] <- "intergenic"
rownames(dat_kwtest) <- NULL

### EXTRACT BY FEATURES ---
dm_kwtest_prom <- subset(dat_kwtest, dat_kwtest$Feature %in% "promoter")
dm_kwtest_intron <- subset(dat_kwtest, dat_kwtest$Feature %in% "intron")
dm_kwtest_intergenic <- subset(dat_kwtest, dat_kwtest$Feature %in% "intergenic")


### GET GRANGES ---
gr_prom <- getGRanges(dat=dm_kwtest_prom)
gr_intron <- getGRanges(dat=dm_kwtest_intron)
gr_intergenic <- getGRanges(dat=dm_kwtest_intergenic)



#############################################################################################################################
### SET SEED ---
set.seed(12345)

### GREAT ENRICHMENT ---
obj_prom <- rGREAT::great(gr=gr_prom, gene_sets="MSigDB:H", tss_source="txdb:hg38", 
                            biomart_dataset = NULL, min_gene_set_size = 5, 
                            mode = "basalPlusExt", basal_upstream = 5000, basal_downstream = 1000, extension = 1000000, 
                            extended_tss = NULL, background = NULL, exclude = "gap", cores = 50, verbose = TRUE)

obj_intron <- rGREAT::great(gr=gr_intron, gene_sets="MSigDB:H", tss_source="txdb:hg38", 
                            biomart_dataset = NULL, min_gene_set_size = 5, 
                            mode = "basalPlusExt", basal_upstream = 5000, basal_downstream = 1000, extension = 1000000, 
                            extended_tss = NULL, background = NULL, exclude = "gap", cores = 50, verbose = TRUE)

obj_intergenic <- rGREAT::great(gr=gr_intergenic, gene_sets="MSigDB:H", tss_source="txdb:hg38", 
                            biomart_dataset = NULL, min_gene_set_size = 5, 
                            mode = "basalPlusExt", basal_upstream = 5000, basal_downstream = 1000, extension = 1000000, 
                            extended_tss = NULL, background = NULL, exclude = "gap", cores = 50, verbose = TRUE)

### PARSE ENRICHMENT DATA ---
dat_prom <- obj_prom@table
dat_intron <- obj_intron@table
dat_intergenic <- obj_intergenic@table

### PROCESS ENRICHMENT DATA ---
df_prom <- processEnrichment(dat=dat_prom, pval=0.05)
df_intron <- processEnrichment(dat=dat_intron, pval=0.05)
df_intergenic <- processEnrichment(dat=dat_intergenic, pval=0.05)




#############################################################################################################################
### GET PLOT  ---
p1 <- plotBar(dat=df_prom, analysis_id="promoter")
p2 <- plotBar(dat=df_intron, analysis_id="intron")
p3 <- plotBar(dat=df_intergenic, analysis_id="intergenic")

### GENERATE PLOT ---
file.plot <- file.path(dir.reproduce_fig, "figure_3_d_f.pdf")
pdf(file.plot, width=3, height=4)
    grid.arrange(p1, p2, p3, ncol=1, nrow=3, heights=c(0.3, 0.4, 0.3))
dev.off()
