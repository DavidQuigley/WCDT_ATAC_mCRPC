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

### DEFINE FILES ----
file.rds_atacseq_masterdata <- file.path(dir.reproduce_data, "wcdt_atacseq_masterdata.rds")
file.rds_rnaseq <- file.path(dir.reproduce_data, "wcdt_rnaseq_atacseq_gene_expression.rds")

### COLOR ---
color_subtypes <- c("#e31a1c","#fb9a99","#ff7f00","#33a02c","#1f78b4")


######################################################################
### LOAD ATAC-SEQ MASTER DATA ---
list.rds_atacseq_masterdata <- readRDS(file=file.rds_atacseq_masterdata)
metadata <- list.rds_atacseq_masterdata$metadata_mcrpc


########## PREPARE RNA-SEQ EXPRESSION DATA -----------------------------------------------------------------------------------
### LOAD RNA-SEQ GENE EXPRESSION OBJECT: TPM ---
cat(format(Sys.time(), "%b %d %X"), "LOADAING R DATA-OBJECT: RNA-seq EXPRESSION ... ", "\n", sep=" ")
    list.rnaseq <- readRDS(file=file.rds_rnaseq)

    # ANNOTATION ---
    gannot <- list.rnaseq$gene_annotation
    gannot <- subset(gannot, gannot$GeneType == "protein_coding")
    gannot <- gannot[-which(stringr::str_detect(gannot$EnsemblID, "PAR_Y")),]

    # EXPRESSION ---
    expr <- as.data.frame(list.rnaseq$gene_expression)
    expr <- subset(expr, expr$FEATURE_ID %in% gannot$EnsemblID)
    rownames(expr) <- gannot$Gene
    expr$FEATURE_ID <- NULL
    colnames(expr) <- stringr::str_replace(colnames(expr), "RP", "BL")

cat(format(Sys.time(), "%b %d %X"), "DONE! ", "\n", sep=" ")

### RE-ARRANGE EXPRESSION DATA ---
expr <- expr[order(rownames(expr), decreasing=FALSE),]
expr <- subset(expr, select=metadata$Sample_ID)

### LOG2 TRANSFORM [log2(TPM + 1)]---
expr <- log2(expr + 1)

########## -------------------------------------------------------------------------------------------------------------------


### SUBSET EXPRESSION DATA BY GENES ---
gexpr <- subset(expr, rownames(expr) == "ZNF263")

### RESHAPE DATA ---
df <- reshape2::melt(t(gexpr))
colnames(df) <- c("SampleID","Gene","Value")



# ADD SUBTYPE AR ---
df$Group <- NA
df$Group[which(df$SampleID %in% metadata$Sample_ID[which(metadata$mCRPC_Subtype == "AR+_NE-")])] <- "AR+_NE-"
df$Group[which(df$SampleID %in% metadata$Sample_ID[which(metadata$mCRPC_Subtype == "ARL_NE-")])] <- "ARL_NE-"
df$Group[which(df$SampleID %in% metadata$Sample_ID[which(metadata$mCRPC_Subtype == "AR+_NE+")])] <- "AR+_NE+"
df$Group[which(df$SampleID %in% metadata$Sample_ID[which(metadata$mCRPC_Subtype == "AR-_NE+")])] <- "AR-_NE+"
df$Group[which(df$SampleID %in% metadata$Sample_ID[which(metadata$mCRPC_Subtype == "AR-_NE-")])] <- "AR-_NE-"

### FACTORIZE ---
df$Group <- factor(df$Group, levels=c("AR+_NE-","ARL_NE-","AR+_NE+","AR-_NE+","AR-_NE-"))


### PLOT ---
p <- ggplot(df, aes(x=Group, y=Value)) +
            geom_boxplot(aes(fill=Group), color="#000000", alpha=0.6, lwd=0.2, outlier.size=0.1, outlier.color="#bdbdbd", na.rm=TRUE) +
            geom_jitter(width=0.1, cex=0.2, color="#525252", alpha=0.5, na.rm=TRUE) +    
            coord_cartesian(ylim=c(2,6)) +
            scale_y_continuous(breaks=seq(2,6, by=1)) +
            scale_fill_manual(values=color_subtypes) +            
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
		    ylab("log2(TPM + 1)") +            
		    xlab("") + 
            ggtitle("ZNF263") 

### PLOT ---
file.plot <- file.path(dir.reproduce_fig, "supplementary_figure_21_d.pdf")
pdf(file.plot, width=3, height=1.5)
    grid.arrange(p, nrow=1, ncol=1)
dev.off()

