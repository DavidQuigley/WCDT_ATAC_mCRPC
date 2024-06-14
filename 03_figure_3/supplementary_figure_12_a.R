###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### DEFINE LIBRARIES ---
library("stringr")
library("data.table")
library("ggplot2")
library("gridExtra")

### DEFINE PATH ---
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/analysis/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

### DEFINE FILES ----
file.rds_atacseq_read_counts_mcrpc <- file.path(dir.reproduce_data, "atacseq_read_counts_mcrpc.rds")
file.kwtest <- file.path(dir.reproduce_data, "wcdt_atacseq_mcrpc_subtypes_KruskalWallisTest_summary.tsv.gz")

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

### FORMAT DIST TO TSS ---
dat_kwtest$distanceToTSS_kb <- abs(dat_kwtest$distanceToTSS/1000)
dat_kwtest$distanceToTSS_kb_log <- log2(dat_kwtest$distanceToTSS_kb + 1)

### ADD LABEL FOR PLOTTING ----
dat_kwtest$Phenotype <- "x"





#############################################################################################################################
# PLOT ---
p <- ggplot(dat_kwtest, aes(x=Phenotype, y=distanceToTSS_kb_log)) +
            #geom_boxplot(fill="#fdbf6f", color="#000000", alpha=0.6, lwd=0.2, outlier.size=0.1, outlier.color="#bdbdbd", na.rm=TRUE) +
            geom_violin(fill="#fdbf6f", color="#000000", alpha=0.6, lwd=0.2, draw_quantiles=NULL, trim=TRUE) +
            geom_boxplot(fill="#FFFFFF", color="#000000", alpha=1, lwd=0.2, width=0.1, outlier.size=0.1, outlier.color="#bdbdbd", na.rm=TRUE) +
            #geom_jitter(width=0.1, cex=0.2, color="#525252", alpha=0.5, na.rm=TRUE) +    
            coord_cartesian(ylim=c(0,10)) +
            scale_y_continuous(breaks=seq(0,10, by=2), labels=c(0,3,15,63,255,1023)) +
            geom_hline(yintercept=2, color="gray70", size=0.4, linetype=2) +
		    theme(
			    axis.text.x = element_text(size = 5, color="#000000"),
			    axis.text.y = element_text(size = 5, color="#000000"),
			    axis.title = element_text(size = 5, color="#000000"),
			    plot.title = element_text(size = 5, color="#000000", hjust=0.5),
			    panel.grid.major = element_blank(),
			    panel.grid.minor = element_blank(),
			    axis.ticks = element_line(size = 0.2, color="#000000"),	
			    panel.background = element_rect(fill = "#FFFFFF", color = "#000000"),
			    legend.text = element_text(size = 5, color="#000000"),
			    legend.title = element_blank(),
			    legend.key.size = unit(0.2, "cm"),			
			    legend.position="none",
                legend.box = "horizontal") + 
            guides(fill = guide_legend(nrow = 2)) +                   
		    xlab("") +            
		    ylab("distance To TSS (kb)") + 
            #ylab("log2(distToTSS + 1)") + 
            ggtitle("") 

### GENERATE PLOT ---
file.plot <- file.path(dir.reproduce_fig, "supplementary_figure_12_a.pdf")
pdf(file.plot, width=1.5, height=2.8)
    grid.arrange(p, nrow=1, ncol=1)
dev.off()


