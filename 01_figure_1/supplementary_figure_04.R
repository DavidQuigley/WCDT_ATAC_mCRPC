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
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

### DEFINE FILES ----
file.rds_atacseq_masterdata <- file.path(dir.reproduce_data, "wcdt_atacseq_masterdata.rds")

### LOAD ATAC-SEQ MASTER DATA ---
list.rds_atacseq_masterdata <- readRDS(file=file.rds_atacseq_masterdata)

### GET STATS DATA ---
dat_pome <- list.rds_atacseq_masterdata$atacseq_stats_pomerantz2020
dat_pome$Dataset <- ifelse(dat_pome$Phenotype == "Benign", "POME_N", "POME_T")

dat_tcga <- list.rds_atacseq_masterdata$atacseq_stats_tcgaprad
dat_tcga$Phenotype <- "PCa"
dat_tcga$Dataset <- "TCGA"

dat_wcdt <- list.rds_atacseq_masterdata$metadata_mcrpc
dat_wcdt <- subset(dat_wcdt, select=c("Sample_ID","ATACseq_Mapped_Reads","ATACseq_Peak_Counts","ATACseq_FRiP"))
colnames(dat_wcdt) <- c("SampleID","ReadCount","PeakCount","FRiP")
dat_wcdt$Phenotype <- "mCRPC"
dat_wcdt$Dataset <- "WCDT"

### MERGET DATA ---
dm <- rbind(dat_wcdt, dat_tcga, dat_pome)

### FACTORIZE ---
dm$Dataset <- factor(dm$Dataset, levels=rev(c("WCDT","TCGA","POME_T","POME_N")))

### COLOR ---
cpalette <- rev(c("#E41A1C","#377EB8","#377EB8","#ffff99"))



### PLOT ---
p1 <- ggplot(dm, aes(x=ReadCount, y=Dataset)) +
            geom_boxplot(aes(fill=Dataset), color="#000000", alpha=0.8, lwd=0.2, outlier.size=0.1, outlier.color="#bdbdbd", na.rm=TRUE) +
            geom_jitter(width=0.01, cex=0.2, color="#525252", alpha=0.5, na.rm=TRUE) +    
            #coord_cartesian(xlim=c(0,50), ylim=c(0,100)) +
            scale_fill_manual(values=cpalette) +            
		    theme(
			    axis.text.x = element_text(size = 5, color="#000000"),
			    axis.text.y = element_text(size = 5, color="#000000"),
			    axis.title = element_text(size = 6, color="#000000"),
			    plot.title = element_text(size = 8, color="#000000", hjust=0.5),
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
		    xlab("No. of Mapped Reads") + 
            ggtitle("")

### PLOT ---
p2 <- ggplot(dm, aes(x=PeakCount, y=Dataset)) +
            geom_boxplot(aes(fill=Dataset), color="#000000", alpha=0.8, lwd=0.2, outlier.size=0.1, outlier.color="#bdbdbd", na.rm=TRUE) +
            geom_jitter(width=0.01, cex=0.2, color="#525252", alpha=0.5, na.rm=TRUE) +    
            #coord_cartesian(xlim=c(0,50), ylim=c(0,100)) +
            scale_fill_manual(values=cpalette) +            
		    theme(
			    axis.text.x = element_text(size = 5, color="#000000"),
			    axis.text.y = element_text(size = 5, color="#000000"),
			    axis.title = element_text(size = 6, color="#000000"),
			    plot.title = element_text(size = 8, color="#000000", hjust=0.5),
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
		    xlab("No. of ATAC-seq Peaks") + 
            ggtitle("")
             
### PLOT ---
p3 <- ggplot(dm, aes(x=FRiP, y=Dataset)) +
            geom_boxplot(aes(fill=Dataset), color="#000000", alpha=0.8, lwd=0.2, outlier.size=0.1, outlier.color="#bdbdbd", na.rm=TRUE) +
            geom_jitter(width=0.01, cex=0.2, color="#525252", alpha=0.5, na.rm=TRUE) +    
            #coord_cartesian(xlim=c(0,50), ylim=c(0,100)) +
            scale_fill_manual(values=cpalette) +            
		    theme(
			    axis.text.x = element_text(size = 5, color="#000000"),
			    axis.text.y = element_text(size = 5, color="#000000"),
			    axis.title = element_text(size = 6, color="#000000"),
			    plot.title = element_text(size = 8, color="#000000", hjust=0.5),
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
		    xlab("FRiP") + 
            ggtitle("")

### PLOT ---
file.plot <- file.path(dir.reproduce_fig, "supplementary_figure_04.pdf")
pdf(file.plot, width=6, height=2)
    grid.arrange(p1, p2, p3, nrow=1, ncol=3)
dev.off()


