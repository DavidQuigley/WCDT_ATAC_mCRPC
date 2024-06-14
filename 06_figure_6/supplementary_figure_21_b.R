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

### MKSCC TAYLOR ET AL. 2010 DATA DOWNLOADED FROM CBIOPORTAL ##############
# url: https://cbioportal-datahub.s3.amazonaws.com/prad_mskcc.tar.gz
# DOWNLOAD THE DATA TO dir.taylor2010

### DEFINE FILES ---
file.des <- file.path(dir.taylor2010, "prad_mskcc/data_clinical_sample.txt")
file.expr <- file.path(dir.taylor2010, "prad_mskcc/data_mrna_agilent_microarray.txt")



### FUNCTION: get.normalizeQuantile() ---
get.normalizeQuantile <- function(dat){
	dQNorm <- data.frame(aroma.light::normalizeQuantile(as.matrix(dat)))
    colnames(dQNorm) <- colnames(dat)
	return(dQNorm)
}	



### LOAD CLINICAL DATA ---
des <-  data.table::fread(file=file.des, sep="\t", header=TRUE, nThread=1, skip=4, data.table=FALSE, verbose=FALSE)
des <- subset(des, des$SAMPLE_CLASS == "Tumor")


### LOAD EXPRESSION DATA ---
dat <-  data.table::fread(file=file.expr, sep="\t", header=TRUE, nThread=1, data.table=FALSE, verbose=FALSE)

### GET GENEID INDEX ---
dm_geneid <- data.frame(Entrez_Gene_Id=dat$Entrez_Gene_Id, Index=c(1:nrow(dat)))

### GET EXPRESSION MATRIX ---
expr <- dat
expr$Entrez_Gene_Id <- NULL

ids <- intersect(des$SAMPLE_ID, colnames(expr))
des <- subset(des, des$SAMPLE_ID %in% ids)
des <- des[match(ids, des$SAMPLE_ID),]
expr <- subset(expr, select=ids)

# HGNC symbol	NCBI gene (formerly Entrezgene) ID
# ZNF263	10127


### QUANTILE NORMALIZE ---
expr_dQnorm <- get.normalizeQuantile(dat=expr)



### RESHAPE DATA ---
df <- reshape2::melt(as.matrix(t(expr_dQnorm)))
colnames(df) <- c("SampleID","Gene","Expr")

### EXTRACT EXPRESSION DATA FOR ZNF263 ---
df <- subset(df, df$Gene == dm_geneid$Index[which(dm_geneid$Entrez_Gene_Id == 10127)])
rownames(df) <- NULL
df$sample_type <- des$SAMPLE_TYPE


### FACTORIZE ---
df$sample_type <- factor(df$sample_type, levels=c("Primary","Metastasis"))


### PLOT ---
p1 <- ggplot(df, aes(x=sample_type, y=Expr)) + 
		    geom_boxplot(aes(fill=sample_type), color="#000000", alpha=0.9, lwd=0.3, outlier.size=0.4, outlier.color="#bdbdbd", na.rm=TRUE) +
            geom_jitter(width=0.1, cex=0.4, color="#525252", alpha=0.5, na.rm=TRUE) +     
            scale_fill_manual(values=c("#377EB8","#E41A1C")) +    
            coord_cartesian(ylim=c(7,8.5)) +
			scale_y_continuous(breaks=seq(7,8.5, by=0.25), labels=seq(7,8.5, by=0.25)) +        
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
		    ylab("mRNA expression (log2)") +            
		    xlab("") + 
            ggtitle("") 


### GENERATE PLOT ---
file.plot <- file.path(dir.reproduce_fig, "supplementary_figure_21_b.pdf")
pdf(file.plot, width=1.5, height=2)
    grid.arrange(p1, nrow=1, ncol=1)  
dev.off()

### WILCOXN TEST ---
pvalue_1_2 <- format(suppressWarnings(wilcox.test(df$Expr[which(df$sample_type == "Primary")], df$Expr[which(df$sample_type == "Metastasis")], paired=FALSE, alternative="two.sided", na.action=na.omit))$p.value, scientific=TRUE)

#> pvalue_1_2
#[1] "2.284502e-05"

