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
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)


### DEFINE FILES ---
file.des <- file.path("/data1/projects/WCDT_deepRNAseq/summary_files/202112_deepRNAseq_samples_key_annotations.txt")
file.expr <- file.path("/data1/projects/WCDT_deepRNAseq/summary_files/gene_counts/featurecounts_gene_tpm.tsv")

# ZNF263 "ENSG00000006194.10"

### LOAD DESIGN TABLE DATA ---
des <-  data.table::fread(file=file.des, sep="\t", header=TRUE, nThread=1, data.table=FALSE, verbose=FALSE)

### ADD SAMPLE TYPE ---
des$sample_type <- NA
des$sample_type[which(des$dataset == "WCDT")] <- "mCRPC"
des$sample_type[which(des$dataset == "CPCG")] <- "PCa"
des$sample_type[which( (des$dataset == "PAIR") & (des$disease_type == "localized-PAIR")  )] <- "PCa"
des$sample_type[which( (des$dataset == "PAIR") & (des$disease_type == "normal-PAIR")  )] <- "normal"


### LOAD EXPRESSION DATA ---
sampleids <-  data.table::fread(file=file.expr, sep="\t", header=FALSE, nrow=1, nThread=50, data.table=FALSE, verbose=FALSE)
expr <-  data.table::fread(file=file.expr, sep="\t", header=FALSE, nThread=50, data.table=FALSE, verbose=FALSE)
colnames(expr) <- c("EnsembleID", as.character(sampleids))
rownames(expr) <- expr$EnsembleID
expr$EnsembleID <- NULL

expr <- subset(expr, select=intersect(des$sample_id, colnames(expr)) )


### LOG2 TRANSFORM [log2(TPM + 1)]---
expr <- log2(expr + 1)





### RESHAPE DATA ---
df <- reshape2::melt(as.matrix(t(expr)))
colnames(df) <- c("SampleID","Gene","Expr")

### EXTRACT EXPRESSION DATA FOR ZNF263 ---
df <- subset(df, df$Gene == "ENSG00000006194.10")
rownames(df) <- NULL
df$sample_type <- des$sample_type
df$disease_type <- des$disease_type


### FACTORIZE ---
df$sample_type <- factor(df$sample_type, levels=c("normal","PCa","mCRPC"))
df$disease_type <- factor(df$disease_type, levels=c("normal-PAIR","localized-PAIR","localized-CPCG","adeno","tSCNC"))




### PLOT ---
p1 <- ggplot(df, aes(x=sample_type, y=Expr)) + 
		    geom_boxplot(aes(fill=sample_type), color="#000000", alpha=0.9, lwd=0.3, outlier.size=0.4, outlier.color="#bdbdbd", na.rm=TRUE) +
            geom_jitter(width=0.1, cex=0.4, color="#525252", alpha=0.5, na.rm=TRUE) +     
            scale_fill_manual(values=c("#ffff99","#377EB8","#E41A1C")) +            
		    theme(
			    axis.text.x = element_text(size = 6, color="#000000"),
			    axis.text.y = element_text(size = 6, color="#000000"),
			    axis.title = element_text(size = 6, color="#000000"),
			    plot.title = element_text(size = 8, color="#000000", hjust=0.5),
			    panel.grid.major = element_blank(),
			    panel.grid.minor = element_blank(),
			    axis.ticks = element_line(size = 0.2, color="#000000"),	
			    panel.background = element_rect(fill = "#FFFFFF", color = "#000000"),
			    legend.text = element_text(size = 6, color="#000000"),
			    legend.title = element_blank(),
			    legend.key.size = unit(0.2, "cm"),			
			    legend.position="none",
                legend.box = "horizontal") + 
            guides(fill = guide_legend(nrow = 1)) +                   
		    ylab("mRNA expression ( log2(TPM+1) )") +            
		    xlab("") + 
            ggtitle("") 


### GENERATE PLOT ---
file.plot <- file.path(dir.reproduce_fig, "supplementary_figure_21_c.pdf")
pdf(file.plot, width=2, height=2)
    grid.arrange(p1, nrow=1, ncol=1)  
dev.off()



### WILCOXN TEST ---
pvalue_1_2 <- format(suppressWarnings(wilcox.test(df$Expr[which(df$sample_type == "normal")], df$Expr[which(df$sample_type == "PCa")], paired=FALSE, alternative="two.sided", na.action=na.omit))$p.value, scientific=TRUE)
pvalue_1_3 <- format(suppressWarnings(wilcox.test(df$Expr[which(df$sample_type == "normal")], df$Expr[which(df$sample_type == "mCRPC")], paired=FALSE, alternative="two.sided", na.action=na.omit))$p.value, scientific=TRUE)
pvalue_2_3 <- format(suppressWarnings(wilcox.test(df$Expr[which(df$sample_type == "PCa")], df$Expr[which(df$sample_type == "mCRPC")], paired=FALSE, alternative="two.sided", na.action=na.omit))$p.value, scientific=TRUE)



#> pvalue_1_2
#[1] "3.634679e-01" -----> 0.3
#> pvalue_1_3
#[1] "4.031263e-02" -----> 0.04
#> pvalue_2_3
#[1] "1.643079e-03" -----> 0.001


