###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### LOAD LIBRARIES ---
library("stringr")
library("data.table")
library("reshape2")
library("ggplot2")
library("gridExtra")

### DEFINE PATH ---
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/analysis/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

dir.tfpeaks2gene <- file.path(dir.reproduce_data, "tfpeaks2gene")

### DEFINE FILES ---
file.tf2genes <- file.path(dir.tfpeaks2gene, "wcdt_tf2genes_r_0p4_pval_0p05_distTSS_500kb.tsv.gz")
file.rds_rnaseq <- file.path(dir.reproduce_data, "wcdt_rnaseq_atacseq_gene_expression.rds")
file.rds_genelist <- file.path(dir.reproduce_data, "genelist.rds")
file.utility_functions <- file.path(dir.reproduce_scripts, "00_prepare_data/00_01_utility_functions.R")

### LOAD FUNCTIONS ---
source(file.utility_functions)

### MOTIFS ---
motif_prefix_znf <- "MA0528.1_ZNF263"



############################################################################################################################################
### LOAD GENELIST ---
list.rds_genelist <- readRDS(file=file.rds_genelist)

### LOAD RNA-SEQ DATA ---
list.rnaseq <- readRDS(file=file.rds_rnaseq)
gannot <- list.rnaseq$gene_annotation

### LOAD TF2GENES DATA ---
dat_tf2genes <- data.table::fread(file=file.tf2genes, sep="\t", header=TRUE, nThread=50, data.table=FALSE, verbose=FALSE)

### GET MOTIF PEAKS ---
dat_znf <- dat_tf2genes[which(dat_tf2genes$motif_prefix == motif_prefix_znf),]




############################################################################################################################################
############################################ OVER REPRESNTATION ANALYSIS ###################################################################
### DEFINE GENES ---
list.genesets <- list.rds_genelist$msigdb_pathways$hallmark
genes.queryset <- unique(dat_znf$Gene)
genes.refset <- gannot$Gene  #BACKGROUND GENESET: ALL GENES IN RNA-SEQ DATA 


# GET OVER-REPRESENTATION ANALYSIS (HYPER-GEOMERTIC TEST) ---
dat_ora <- get_over_representation_analysis(list.genesets, genes.queryset, genes.refset, p.threshold=0.0001)



############################################################################################################################################
############################################### PLOT ENRICHMENT ############################################################################
### FUNCTION: getPathwaysData() ---
getPathwaysData <- function(dat){
    # EXTRACT DATA ---
    if(nrow(dat) >= 20){
        dat <- dat[1:20,]
    }

    # ADD PVALUES ---
    pval <- min(setdiff(dat$pvalue, 0))
    index0 <- which(dat$pvalue == 0)

    if(length(index0) != 0){
        dat$pvalue[which(dat$pvalue == 0)] <- pval
    }

    # ADD VALUES ---
    dat$nlogp <- -log10(dat$fdr)
    dat <- dat[order(dat$nlogp, decreasing=FALSE),]
    dat$Category <- factor(dat$Category, levels=dat$Category)

    return(dat)
}

### FUNCTION: plotBar() ----
plotBar <- function(df){
    # PLOT ---
    p <- ggplot(df, aes(x=nlogp, y=Category)) +
            geom_bar(stat="identity", fill="#006094") +
            coord_cartesian(xlim=c(100,250)) +
            scale_x_continuous(breaks=seq(100,250,by=50)) +
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
            xlab("-log10(qvalue)") +
            ylab("") + 
            ggtitle("") 
    
    return(p)
}



### LOAD ENRICHMENT DATA ---
dat_ora_parsed <- getPathwaysData(dat=dat_ora)

### GET PLOT ---
p_bar <- plotBar(df=dat_ora_parsed)

### PLOT ---
file.plot <- file.path(dir.reproduce_fig, "figure_6_b.pdf")
pdf(file.plot, width=2.8, height=2)
    grid.arrange(p_bar, nrow=1, ncol=1)
dev.off()



