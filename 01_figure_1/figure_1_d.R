###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### DEFINE LIBRARIES ---
library("stringr")
library("data.table")
library("dplyr")
library("reshape2")
library("ggplot2")
library("gridExtra")

### DEFINE PATH ---
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/analysis/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

### DEFINE FILES ----
file.rds_atacseq_chromatin_variants <- file.path(dir.reproduce_data, "atacseq_chromatin_variants_benign_localizedpca_mcrpc_adeno_nepc.rds")
file.rds_atacseq_read_counts_combined <- file.path(dir.reproduce_data, "atacseq_read_counts_combined.rds")


### LOAD ATAC-SEQ READ COUNTS DATA ---
list.rds_atacseq_read_counts_combined <- readRDS(file=file.rds_atacseq_read_counts_combined)

### GET FEATURE ANNOTATION DATA ---
dat.annot <- list.rds_atacseq_read_counts_combined$feature_annotation
dat.annot <- subset(dat.annot, select=c("FeatureID","Feature"))


### LOAD ATAC-SEQ CHROMATIN VARIANTS DATA ---
list.rds_atacseq_chromatin_variants <- readRDS(file=file.rds_atacseq_chromatin_variants)

### GET CLUSTERED FEATURES ---
dat <- list.rds_atacseq_chromatin_variants$atacseq_diffpeaks_cluster_features
dat <- subset(dat, dat$ClusterName == "mcrpc_adeno")



### MERGE DATA ---
dat <- merge(dat, dat.annot, by="FeatureID")
dat <- subset(dat, dat$Feature %in% c("Promoter","Intron","Distal Intergenic"))
dat$Feature[which(dat$Feature == "Distal Intergenic")] <- "Intergenic"


### GET PEAK FEATURe COUNT ---
dm <- dat %>%
        dplyr::select(Feature, ClusterName, FeatureID) %>%
        dplyr::group_by(Feature, ClusterName) %>%
        dplyr::count(Feature, ClusterName, sort=FALSE, name="Freq")

dm$Percentage <- round( (dm$Freq/sum(dm$Freq)) * 100, 2)


### FACTORIZE ---
dm$Feature <- factor(dm$Feature, levels=c("Promoter","Intron","Intergenic"))

### PLOT ---
cpallete <- c("#6a3d9a","#b15928","#e31a1c")
p <- ggplot(dm, aes(x=Feature, y=Percentage, fill=Feature)) + 
		geom_bar(stat="identity", color="#000000", width=0.8, size=0.2) + 
		scale_fill_manual(values=cpallete) +
		coord_cartesian(ylim=c(0,50)) +
		scale_y_continuous(breaks=seq(0,50, by=10)) + 
		theme(
			axis.text.x = element_text(size = 6, color="#000000"),
			axis.text.y = element_text(size = 6, color="#000000"),
			axis.title = element_text(size = 8, color="#000000"),
			plot.title = element_text(size = 10, color="#000000", hjust=0),
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			axis.ticks = element_line(size = 0.2, colour="#000000"),	
			panel.background = element_rect(fill = "#FFFFFF", colour = "#000000"),           
			legend.text = element_text(size = 5, color="#000000"),
			legend.title = element_blank(),
			legend.key.size = unit(0.3, "cm"),			
			legend.position="none") + 
			ylab("% of chromatin variants") +
			xlab("") + 
            ggtitle("") 

### GENERATE PLOT ---
file.plot <- file.path(dir.reproduce_fig, "figure_1_d.pdf")
pdf(file.plot, height=2, width=2)
    grid.arrange(p, nrow=1, ncol=1)
dev.off()


