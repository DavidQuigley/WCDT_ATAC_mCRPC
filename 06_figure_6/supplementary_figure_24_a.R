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
dat_znf <- dat_znf[!duplicated(dat_znf),]
dat_znf <- subset(dat_znf, dat_znf$Corr > 0)
rownames(dat_znf) <- NULL





############################################################################################################################################
############################################ OVER REPRESNTATION ANALYSIS ###################################################################
### DEFINE GENES ---
list.genesets <- list.rds_genelist$msigdb_pathways$hallmark
genes.queryset <- unique(dat_znf$Gene)
genes.refset <- gannot$Gene  #BACKGROUND GENESET: ALL GENES IN RNA-SEQ DATA 


# GET OVER-REPRESENTATION ANALYSIS (HYPER-GEOMERTIC TEST) ---
dat_ora <- get_over_representation_analysis(list.genesets, genes.queryset, genes.refset, p.threshold=0.0001)
dat_enrich <- dat_ora[1:20,]






############################################################################################################################################
### EXTRACT GENES CORR BY PATHWAY ---
list.dm <- list()
for(i in 1:nrow(dat_enrich)){
    item <- dat_enrich$Category[i]
    genes_item <- stringr::str_split( dat_enrich$overlap.genes[which(dat_enrich$Category == item)], ":")[[1]]
    
    d <- subset(dat_znf, dat_znf$Gene %in% genes_item)
    d$Analysis <- item

    list.dm[[item]] <- d
}

### AGGREAGATE DATA ---
dm <- do.call(rbind.data.frame, list.dm)
rownames(dm) <- NULL





############################################################################################################################################
### FACTORIZE ---
dm$Analysis <- factor(dm$Analysis, levels=rev(dat_enrich$Category))

### BOX-PLOT ---
p_box <- ggplot(dm, aes(x=Corr, y=Analysis)) +
            geom_boxplot(fill="#e31a1c", color="#000000", alpha=0.6, lwd=0.15, outlier.size=0.1, outlier.color="#bdbdbd", na.rm=TRUE) +
            #geom_jitter(width=0.1, cex=0.2, color="#525252", alpha=0.5, na.rm=TRUE) +    
            coord_cartesian(xlim=c(-0.6,0.8)) +
            scale_x_continuous(breaks=round(seq(-0.6, 0.8, by=0.2), 1), labels=round(seq(-0.6, 0.8, by=0.2), 1) ) +
            #scale_fill_manual(values=cpalette) +            
		    theme(
			    axis.text.x = element_text(size = 5, color="#000000"),
			    axis.text.y = element_text(size = 5, color="#000000"),
			    axis.title = element_text(size = 5, color="#000000"),
			    plot.title = element_text(size = 5, color="#000000", hjust=0.5),
			    panel.grid.major = element_blank(),
			    panel.grid.minor = element_blank(),
			    axis.ticks = element_line(size = 0.2, color="#000000"),	
                strip.text.x = element_text(size=5, color="#000000"),
                strip.text.y = element_text(size=5, color="#000000", angle=0, hjust=0),
                strip.background = element_rect(fill="#FFFFFF", color="#FFFFFF"),
                panel.background = element_rect(fill="#FFFFFF", color="#000000"),
			    legend.text = element_text(size = 5, color="#000000"),
			    legend.title = element_blank(),
			    legend.key.size = unit(0.2, "cm"),			
			    legend.position="none",
                legend.box = "horizontal") + 
            guides(fill = guide_legend(nrow = 2)) +                   
		    ylab("") +            
		    xlab("Spearman's Correlation Coefficient (R)") + 
            ggtitle("") 

### GENERATE PLOT ---
file.plot <- file.path(dir.reproduce_fig, "supplementary_figure_24_a.pdf")
pdf(file.plot, width=3, height=2)
    grid.arrange(p_box, nrow=1, ncol=1)   
dev.off()
