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


######################################################################
### LOAD ATAC-SEQ MASTER DATA ---
list.rds_atacseq_masterdata <- readRDS(file=file.rds_atacseq_masterdata)

### GET METADATA ---
metadata <- list.rds_atacseq_masterdata$metadata_mcrpc
des <- subset(metadata, select=c("Sample_ID","mCRPC_Subtype","Biopsy_Site","ATACseq_Peak_Counts","ATACseq_FRiP"))
colnames(des) <- c("Sample_ID","Subtype","Biopsy_Site","Peaks","FRiP")


### PEAK NUMBERS IN 10000 ---
des$Peaks <- des$Peaks/10000

### FACTORIZE ---
des$Subtype <- factor(des$Subtype, levels=c("AR+_NE-","ARL_NE-","AR+_NE+","AR-_NE+","AR-_NE-"))

### PLOT ---
p1 <- ggplot(des, aes(x=Subtype, y=Peaks)) + 
		    geom_boxplot(aes(fill=Subtype), color="#000000", alpha=0.9, lwd=0.3, outlier.size=0.4, outlier.color="#bdbdbd", na.rm=TRUE) +
            geom_jitter(width=0.1, cex=0.4, color="#525252", alpha=0.5, na.rm=TRUE) +        
            coord_cartesian(ylim=c(2, 16)) +
            scale_y_continuous(breaks=seq(2,16, by=2), labels=seq(2,16, by=2)) +
            scale_fill_manual(values=c("#e31a1c","#fb9a99","#ff7f00","#33a02c","#1f78b4")) +            
		    theme(
			    axis.text.x = element_text(size = 5, color="#000000", angle=45, hjust=0.5, vjust=0.5),
			    axis.text.y = element_text(size = 5, color="#000000"),
			    axis.title = element_text(size = 5, color="#000000"),
			    plot.title = element_text(size = 5, color="#000000", hjust=0.5),
			    panel.grid.major = element_blank(),
			    panel.grid.minor = element_blank(),
			    axis.ticks = element_line(size = 0.2, colour="#000000"),	
			    panel.background = element_rect(fill = "#FFFFFF", colour = "#000000"),
			    legend.text = element_text(size = 5, color="#000000"),
			    legend.title = element_blank(),
			    legend.key.size = unit(0.2, "cm"),			
			    legend.position="none",
                legend.box = "horizontal") + 
            guides(fill = guide_legend(nrow = 1)) +                   
		    ylab(  expression("No. of ATAC-seq Peaks in ( 10"^4*" )") ) +                
		    xlab("mCRPC Subtypes") + 
            ggtitle("") 

### PLOT ---
p2 <- ggplot(des, aes(x=Subtype, y=FRiP)) + 
		    geom_boxplot(aes(fill=Subtype), color="#000000", alpha=0.9, lwd=0.3, outlier.size=0.4, outlier.color="#bdbdbd", na.rm=TRUE) +
            geom_jitter(width=0.1, cex=0.4, color="#525252", alpha=0.5, na.rm=TRUE) +     
            coord_cartesian(ylim=c(0,0.3)) +
            scale_y_continuous(breaks=seq(0,0.3, by=0.05), labels=seq(0,0.3, by=0.05) ) +
            scale_fill_manual(values=c("#e31a1c","#fb9a99","#ff7f00","#33a02c","#1f78b4")) +            
		    theme(
			    axis.text.x = element_text(size = 6, color="#000000", angle=45, hjust=0.5, vjust=0.5),
			    axis.text.y = element_text(size = 6, color="#000000"),
			    axis.title = element_text(size = 6, color="#000000"),
			    plot.title = element_text(size = 6, color="#000000", hjust=0.5),
			    panel.grid.major = element_blank(),
			    panel.grid.minor = element_blank(),
			    axis.ticks = element_line(size = 0.2, colour="#000000"),	
			    panel.background = element_rect(fill = "#FFFFFF", colour = "#000000"),
			    legend.text = element_text(size = 6, color="#000000"),
			    legend.title = element_blank(),
			    legend.key.size = unit(0.2, "cm"),			
			    legend.position="none",
                legend.box = "none") + 
            guides(fill = guide_legend(nrow = 1)) +                   
		    ylab("FRiP") +            
		    xlab("mCRPC Subtypes") + 
            ggtitle("") 

### GENERATE PLOT ---
file.plot <- file.path(dir.reproduce_fig, "supplementary_figure_11.pdf")
pdf(file.plot, height=2, width=4)
    grid.arrange(p1, p2, nrow=1, ncol=2)
dev.off()

