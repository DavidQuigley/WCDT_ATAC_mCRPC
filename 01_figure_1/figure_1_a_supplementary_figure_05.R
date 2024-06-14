###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### DEFINE LIBRARIES ---
library("stringr")
library("data.table")
library("scatterplot3d")
library("ggplot2")
library("gridExtra")

### DEFINE PATH ---
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/analysis/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

### DEFINE FILES ----
file.rds_atacseq_masterdata <- file.path(dir.reproduce_data, "wcdt_atacseq_masterdata.rds")
file.rds_atacseq_read_counts_combined <- file.path(dir.reproduce_data, "atacseq_read_counts_combined.rds")
file.utility_functions <- file.path(dir.reproduce_scripts, "00_prepare_data/00_01_utility_functions.R")

### LOAD FUNCTIONS ---
source(file.utility_functions)

### LOAD ATAC-SEQ MASTER DATA ---
list.rds_atacseq_masterdata <- readRDS(file=file.rds_atacseq_masterdata)

### LOAD ATAC-SEQ READ COUNTS DATA ---
list.rds_atacseq_read_counts_combined <- readRDS(file=file.rds_atacseq_read_counts_combined)
dat_read_counts_norm <- list.rds_atacseq_read_counts_combined$feature_counts_norm
rownames(dat_read_counts_norm) <- dat_read_counts_norm$feature_id
dat_read_counts_norm$feature_id <- NULL


### FUNCTION: getPCAdata() ---
getPCAdata <- function(mat, topVarPercentile){
    # EXTRACT FEATURES --- # LOADED FROM FILE.utility_functions #
    list.features_filter <- featurefilter(mydata=mat, percentile=topVarPercentile)
    mat <- list.features_filter$filtered_data
    mat <- t(mat)

    # PERFORM A PCA ON THE DATA FOR THE SELECTED FEATURES
    pca <- prcomp(mat)

    # THE CONTRIBUTION TO THE TOTAL VARIANCE FOR EACH COMPONENT
    percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
    percentVar <- percentVar[1:3]
    percentVar <- round(percentVar * 100, 1)
    names(percentVar) <- paste("PC", c(1:3), sep="")

    # PREPARE DATA ---
    d <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], PC3=pca$x[,3])
    
    # COMPILE DATA ---
    list.output <- list(d, percentVar)

    return(list.output)
}

### SET SEED ---
set.seed(12345)

### GET 20% MOST VARIABLE FEATURES ----
list.pca <- getPCAdata(mat=dat_read_counts_norm, topVarPercentile=20)

### EXTRACT DATA FROM PCA ---
df <- list.pca[[1]]
dm <- cbind(df, list.rds_atacseq_masterdata$metadata_combined)



#######################################################################################
### PLOT: PCA 3D FINALIZED ----
file.plot <- file.path(dir.reproduce_fig, "figure_1_a.pdf")
pdf(file.plot, width=3.5, height=3.5)   

    scatterplot3d::scatterplot3d(dm$PC1,    # x axis#
				                dm$PC2,     # y axis
				                dm$PC3,     # z axis
				                bg=dm$Color, color=rep("#000000", nrow(dm)),
				                highlight.3d=FALSE, pch=21, type="p",
				                cex.symbols=0.7, cex.axis=0.3, cex.lab=0.3, 
                                #xlim=c(-300,200), ylim=c(-300,100), zlim=c(-100,300),
				                angle=60, scale.y=1, 
				                grid=TRUE, box=TRUE,
				                main="", xlab="PC1", ylab="PC2", zlab="PC3")

dev.off()




#### VARIANCE CAPTURED ---
#> list.pca[[2]]
# PC1  PC2  PC3
#17.0 13.7  6.8

### VARIANCE DATA ---
dm_var <- data.frame(pc=names(list.pca[[2]]), var=as.numeric(list.pca[[2]]))
dm_var$pc <- factor(dm_var$pc, levels=c("PC1","PC2","PC3"))

### GET PLOT: PC VARIANCE ---
p_var <- ggplot(dm_var, aes(x=pc, y=var)) + 
            geom_bar(stat="identity", color=NA, fill="#000000",  width=0.8, size=0.5) + 
            coord_cartesian(ylim=c(0,18)) +
            scale_y_continuous(breaks=seq(0,18,by=3)) +
            theme(
                axis.text.x = element_text(size = 6, color="#000000", angle=90, hjust=1, vjust=0.5),
                axis.text.y = element_text(size = 6, color="#000000"),
                axis.title = element_text(size = 8, color="#000000"),
                plot.title = element_text(size = 8, color="#000000", hjust=0.5),
                panel.grid.major.y = element_blank(),
                panel.grid.major.x = element_blank(),
                panel.grid.minor.x = element_blank(),
                panel.grid.minor.y = element_blank(),
                axis.ticks = element_line(size=0.2, color="#000000"), 
                strip.text = element_text(size=10, color="#000000"),
                strip.background = element_rect(fill="#FFFFFF", color="#FFFFFF"),
                panel.background = element_rect(fill="#FFFFFF", color="#000000"),
                legend.text = element_text(size = 10, color="#000000"),
                legend.title = element_blank(),
                legend.key.size = unit(0.5, "cm"),
                legend.position = "none") +
            ylab("% Variance in PC") +
            xlab("") + 
            ggtitle("") 


### GENERATE PLOT ---
file.plot <- file.path(dir.reproduce_fig, "supplementary_figure_05_a.pdf")
pdf(file.plot, width=2, height=2)
    grid.arrange(p_var, nrow=1, ncol=1)
dev.off()


#######################################################################################
dm$Phenotype <- factor(dm$Phenotype, levels=c("Benign","PCa","mCRPC-Adeno","mCRPC-SmallCell","mCRPC-SmallCell_Cejas"))
cpalette <- c("#33a02c","#377EB8","#E41A1C","#ffff99","#bdbdbd")

### PLOT: PCA 2D ----
p_pca2d_1 <- ggplot(dm, aes(x=PC1, y=PC2)) + 
                geom_point(aes(fill=Phenotype), color="#000000", shape = 21, stroke=0.5, size=1.25, alpha=1, na.rm=TRUE) +
                scale_fill_manual(values=cpalette) +
                #scale_color_manual(values=cpalette) +
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
                ylab("PCA2") +
                xlab("PCA1") + 
                ggtitle("") 


### PLOT: PCA 2D ----
p_pca2d_2 <- ggplot(dm, aes(x=PC1, y=PC3)) + 
                geom_point(aes(fill=Phenotype), color="#000000", shape = 21, stroke=0.5, size=1.25, alpha=1, na.rm=TRUE) +
                scale_fill_manual(values=cpalette) +
                #scale_color_manual(values=cpalette) +
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
                ylab("PCA3") +
                xlab("PCA1") + 
                ggtitle("") 


### PLOT: PCA 2D ----
p_pca2d_3 <- ggplot(dm, aes(x=PC2, y=PC3)) + 
                geom_point(aes(fill=Phenotype), color="#000000", shape = 21, stroke=0.5, size=1.25, alpha=1, na.rm=TRUE) +
                scale_fill_manual(values=cpalette) +
                #scale_color_manual(values=cpalette) +
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
                ylab("PCA3") +
                xlab("PCA2") + 
                ggtitle("") 

### GENERATE PLOT ---
file.plot <- file.path(dir.reproduce_fig, "supplementary_figure_05_b_d.pdf")
pdf(file.plot, width=6, height=2)
    grid.arrange(p_pca2d_1, p_pca2d_2, p_pca2d_3, nrow=1, ncol=3)
dev.off()
