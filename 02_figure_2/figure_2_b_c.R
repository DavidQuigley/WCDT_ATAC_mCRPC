###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### DEFINE LIBRARIES ---
library("stringr")
library("data.table")
library("singscore")
library("ggplot2")
library("gridExtra")

### DEFINE PATH ---
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

### DEFINE FILES ----
file.rds_atacseq_masterdata <- file.path(dir.reproduce_data, "wcdt_atacseq_masterdata.rds")
file.cluster_membership <- file.path(dir.reproduce_data, "wcdt_atacseq_unsupervised_cluster_membership.tsv") 
file.rds_genelist <- file.path(dir.reproduce_data, "genelist.rds")
file.rds_rnaseq <- file.path(dir.reproduce_data, "wcdt_rnaseq_atacseq_gene_expression.rds")
file.utility_functions <- file.path(dir.reproduce_scripts, "00_prepare_data/00_01_utility_functions.R")

### GENES FROM genes_nepc ---
genes_nepc_up <- c("SVOP","KCNB2","RGS7","KCND2","BRINP1","GNAO1",
					"PCSK1","ST8SIA3","SEZ6","SCG3","NRSN1","TRIM9",
					"LRRC16B","SOGA3","AURKA","DNMT1","EZH2","MYCN",
					"ETV5","GPX2","SYT11","PROX1","ASXL3","JAKMIP2")

genes_nepc_dn <- c("NUP93","FOXP1","HERPUD1","SLC25A37","RAB27B","RGS10",
					"CCND1","RB1","MAPKAPK3","MYH9","MMP2","PIEZO1",
					"PSCA","UPK2","EFNA4","OPHN1","EVPL","GATA2",
					"SLC44A4","HOXB13","EPN3","SPDEF","RIPK2","KLK4",
					"AR","KLK3","NKX3-1")

######################################################################
### LOAD ATAC-SEQ MASTER DATA ---
list.rds_atacseq_masterdata <- readRDS(file=file.rds_atacseq_masterdata)
metadata <- list.rds_atacseq_masterdata$metadata_mcrpc

######################################################################
### LOAD GENELIST DATA ---
list.rds_genelist <- readRDS(file=file.rds_genelist)
genelist.ar <- list.rds_genelist$msigdb_pathways$hallmark$HALLMARK_ANDROGEN_RESPONSE 
genes_nepc <- list.rds_genelist$genes_nepc

### LOAD CLUSTER MEMBERSHIP DATA: FROM UNSUPERVISED CLUSTERING ---
dat.clust <- data.table::fread(file=file.cluster_membership, sep="\t", header=TRUE, nThread=1,  data.table=FALSE, verbose=FALSE)



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





######################################################################
### RANK GENE EXPRESSION DATA ---
gexpr_rank <- singscore::rankGenes(expreMatrix=as.matrix(expr))

### COMPUTE GENESET SCORES ---
dat.score_ar <- singscore::simpleScore(rankData=gexpr_rank, upSet=genelist.ar, knownDirection=TRUE, centerScore=TRUE, subSamples=NULL)
dat.score_ne <- singscore::simpleScore(rankData=gexpr_rank, upSet=genes_nepc_up, downSet=genes_nepc_dn, knownDirection=TRUE, centerScore=TRUE, subSamples=NULL)

### COMBINE DATA ---
dat_score   <- data.frame(SampleID=rownames(dat.score_ar),
                            ARscore=dat.score_ar$TotalScore,
                            NEscore=dat.score_ne$TotalScore)


### SAVE OBJECT TO RDATA FILE ---
file.rds_rnaseq_ar_nepc_scores <- file.path(dir.reproduce_data, "wcdt_rnaseq_ar_nepc_scores.rds")
saveRDS(object=dat_score, file=file.rds_rnaseq_ar_nepc_scores)


### MERGE SCORE AND CLUSTER DATA ---
dm <- merge(dat_score, dat.clust, by="SampleID", sort=FALSE)

### FACTORIZE ---
dm$Cluster <- factor(dm$Cluster, levels=c("clust_1","clust_2","clust_3"))



### PLOT ---
p1 <- ggplot(dm, aes(x=Cluster, y=ARscore)) + 
		    geom_boxplot(aes(fill=Cluster), color="#000000", alpha=0.9, lwd=0.3, outlier.size=0.4, outlier.color="#bdbdbd", na.rm=TRUE) +
            geom_jitter(width=0.1, cex=0.4, color="#525252", alpha=0.5, na.rm=TRUE) +     
            coord_cartesian(ylim=c(0.2,0.45)) +
			scale_y_continuous(breaks=seq(0.2,0.45, by=0.05), labels=seq(0.2,0.45, by=0.05)) +
            scale_fill_manual(values=c("#33a02c","#fb9a99","#e31a1c")) +            
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
		    ylab("AR score") +            
		    xlab("ATAC-seq Clusters") + 
            ggtitle("") 


p2 <- ggplot(dm, aes(x=Cluster, y=NEscore)) + 
		    geom_boxplot(aes(fill=Cluster), color="#000000", alpha=0.9, lwd=0.3, outlier.size=0.4, outlier.color="#bdbdbd", na.rm=TRUE) +
            geom_jitter(width=0.1, cex=0.4, color="#525252", alpha=0.5, na.rm=TRUE) +     
            coord_cartesian(ylim=c(-0.6,0.4)) +
			scale_y_continuous(breaks=round(seq(-0.6,0.4, by=0.2),2), labels=round(seq(-0.6,0.4, by=0.2), 2) ) +
            scale_fill_manual(values=c("#33a02c","#fb9a99","#e31a1c")) +            
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
		    ylab("NE score") +            
		    xlab("ATAC-seq Clusters") + 
            ggtitle("") 

### GENERATE PLOT ---
file.plot <- file.path(dir.reproduce_fig, "figure_2_b_c.pdf")
pdf(file.plot, height=3.5, width=2)
    grid.arrange(p1, p2, nrow=2, ncol=1)
dev.off()




######################################################################

### WILCOXN TEST ---
pvalue_1_2 <- format(suppressWarnings(wilcox.test(dm$ARscore[which(dm$Cluster == "clust_1")], dm$ARscore[which(dm$Cluster == "clust_2")], paired=FALSE, alternative="two.sided", na.action=na.omit))$p.value, scientific=TRUE)
pvalue_1_3 <- format(suppressWarnings(wilcox.test(dm$ARscore[which(dm$Cluster == "clust_1")], dm$ARscore[which(dm$Cluster == "clust_3")], paired=FALSE, alternative="two.sided", na.action=na.omit))$p.value, scientific=TRUE)
pvalue_2_3 <- format(suppressWarnings(wilcox.test(dm$ARscore[which(dm$Cluster == "clust_2")], dm$ARscore[which(dm$Cluster == "clust_3")], paired=FALSE, alternative="two.sided", na.action=na.omit))$p.value, scientific=TRUE)

#> pvalue_1_2
#[1] "7.655229e-04"
#> pvalue_1_3
#[1] "3.065923e-08"
#> pvalue_2_3
#[1] "4.037678e-04"

### WILCOXN TEST ---
pvalue_1_2 <- format(suppressWarnings(wilcox.test(dm$NEscore[which(dm$Cluster == "clust_1")], dm$NEscore[which(dm$Cluster == "clust_2")], paired=FALSE, alternative="two.sided", na.action=na.omit))$p.value, scientific=TRUE)
pvalue_1_3 <- format(suppressWarnings(wilcox.test(dm$NEscore[which(dm$Cluster == "clust_1")], dm$NEscore[which(dm$Cluster == "clust_3")], paired=FALSE, alternative="two.sided", na.action=na.omit))$p.value, scientific=TRUE)
pvalue_2_3 <- format(suppressWarnings(wilcox.test(dm$NEscore[which(dm$Cluster == "clust_2")], dm$NEscore[which(dm$Cluster == "clust_3")], paired=FALSE, alternative="two.sided", na.action=na.omit))$p.value, scientific=TRUE)


#> pvalue_1_2
#[1] "4.989399e-04"
#> pvalue_1_3
#[1] "6.790254e-05"
#> pvalue_2_3
#[1] "5.229901e-01"

