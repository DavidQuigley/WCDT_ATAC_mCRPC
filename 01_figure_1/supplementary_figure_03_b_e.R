###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### DEFINE LIBRARIES ---
library("stringr")
library("data.table")
library("dplyr")
library("ggplot2")
library("gridExtra")

### DEFINE PATH ---
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/analysis/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

### DEFINE FILES ----
file.rds_atacseq_masterdata <- file.path(dir.reproduce_data, "wcdt_atacseq_masterdata.rds")

### FUNCTION: getCorr() ---
getCorr <- function(x, y){
    correlation.test <- cor.test(x, y, method="spearman", na.action = "na.exclude")
	r <- round(as.numeric(correlation.test$estimate),2)
    return(r)
}

### LOAD ATAC-SEQ MASTER DATA ---
list.rds_atacseq_masterdata <- readRDS(file=file.rds_atacseq_masterdata)

### GET METADATA ---
metadata <- subset(list.rds_atacseq_masterdata$metadata_mcrpc, select=c("Sample_ID","ATACseq_Mapped_Reads","ATACseq_Peak_Counts","ATACseq_TSS_Enrich_Score","ATACseq_FRiP"))
metadata_wgs <- subset(list.rds_atacseq_masterdata$stats_wgs_mcrpc, select=c("SampleID","Purity"))

### GET SPEARMAN'S CORRELATION ---
r1 <- getCorr(x=metadata$ATACseq_Peak_Counts, y=metadata$ATACseq_FRiP)
r2 <- getCorr(x=metadata$ATACseq_Peak_Counts, y=metadata$ATACseq_TSS_Enrich_Score)
r3 <- getCorr(x=metadata$ATACseq_FRiP, y=metadata$ATACseq_TSS_Enrich_Score)

cor.test(x=metadata$ATACseq_Peak_Counts, y=metadata$ATACseq_FRiP, method="spearman", na.action = "na.exclude")
cor.test(x=metadata$ATACseq_Peak_Counts, y=metadata$ATACseq_TSS_Enrich_Score, method="spearman", na.action = "na.exclude")
cor.test(x=metadata$ATACseq_FRiP, y=metadata$ATACseq_TSS_Enrich_Score, method="spearman", na.action = "na.exclude")

#> cor.test(x=metadata$ATACseq_Peak_Counts, y=metadata$ATACseq_FRiP, method="spearman", na.action = "na.exclude")
#
#        Spearman's rank correlation rho
#
#data:  metadata$ATACseq_Peak_Counts and metadata$ATACseq_FRiP
#S = 12119, p-value = 5.781e-16
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#      rho
#0.7879642
#
#Warning message:
#In cor.test.default(x = metadata$ATACseq_Peak_Counts, y = metadata$ATACseq_FRiP,  :
#  Cannot compute exact p-value with ties


#> cor.test(x=metadata$ATACseq_Peak_Counts, y=metadata$ATACseq_TSS_Enrich_Score, method="spearman", na.action = "na.exclude")
#
#        Spearman's rank correlation rho
#
#data:  metadata$ATACseq_Peak_Counts and metadata$ATACseq_TSS_Enrich_Score
#S = 44692, p-value = 0.06989
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#      rho
#0.2180562



#> cor.test(x=metadata$ATACseq_FRiP, y=metadata$ATACseq_TSS_Enrich_Score, method="spearman", na.action = "na.exclude")
#
#        Spearman's rank correlation rho
#
#data:  metadata$ATACseq_FRiP and metadata$ATACseq_TSS_Enrich_Score
#S = 19229, p-value = 3.843e-10
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#      rho
#0.6635635
#
##Warning message:
#In cor.test.default(x = metadata$ATACseq_FRiP, y = metadata$ATACseq_TSS_Enrich_Score,  :
#  Cannot compute exact p-value with ties


### PLOT 1 ---
p1 <- ggplot(metadata, aes(x=ATACseq_Peak_Counts, y=ATACseq_FRiP)) + 
            stat_smooth(method="lm", geom = "smooth", formula = y ~ x, position = "identity", size=0.5, fullrange = FALSE, se=TRUE, na.rm=TRUE) +	
            geom_point(fill="#e41a1c", color="#e41a1c", shape = 21, stroke=0.5, size=0.8, alpha=0.7, na.rm=TRUE) +
            coord_cartesian(xlim=c(20000, 160000), ylim=c(0, 0.3)) +
            scale_x_continuous(breaks=seq(20000,160000, by=20000), labels=seq(20000,160000, by=20000)) +
            scale_y_continuous(breaks=seq(0,0.3, by=0.05), labels=seq(0,0.3, by=0.05) ) +
            annotate("text", x = 120000, y = 0.025, label = paste("R = ", r1, sep=""), size=2.5) +	
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
            ylab("FRiP") +
            xlab("ATAC-seq Peak Counts") + 
            ggtitle("") 


### PLOT 2 ---
p2 <- ggplot(metadata, aes(x=ATACseq_Peak_Counts, y=ATACseq_TSS_Enrich_Score)) + 
            stat_smooth(method="lm", geom = "smooth", formula = y ~ x, position = "identity", size=0.5, fullrange = FALSE, se=TRUE, na.rm=TRUE) +	
            geom_point(fill="#e41a1c", color="#e41a1c", shape = 21, stroke=0.5, size=0.8, alpha=0.7, na.rm=TRUE) +
            coord_cartesian(xlim=c(20000, 160000), ylim=c(0, 16)) +
            scale_x_continuous(breaks=seq(20000,160000, by=20000), labels=seq(20000,160000, by=20000)) +
            scale_y_continuous(breaks=seq(0,16, by=4), labels=seq(0,16, by=4) ) +
            annotate("text", x = 120000, y = 1, label = paste("R = ", r2, sep=""), size=2.5) +
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
            ylab("TSS Enrich Score") +
            xlab("ATAC-seq Peak Counts") + 
            ggtitle("") 

### PLOT 3 ---
p3 <- ggplot(metadata, aes(x=ATACseq_FRiP, y=ATACseq_TSS_Enrich_Score)) + 
            stat_smooth(method="lm", geom = "smooth", formula = y ~ x, position = "identity", size=0.5, fullrange = FALSE, se=TRUE, na.rm=TRUE) +	
            geom_point(fill="#e41a1c", color="#e41a1c", shape = 21, stroke=0.5, size=0.8, alpha=0.7, na.rm=TRUE) +
            coord_cartesian(xlim=c(0, 0.3), ylim=c(0, 16)) +
            scale_x_continuous(breaks=seq(0,0.3, by=0.05), labels=seq(0,0.3, by=0.05) ) +
            scale_y_continuous(breaks=seq(0,16, by=4), labels=seq(0,16, by=4) ) +
            annotate("text", x = 0.2, y = 1, label = paste("R = ", r3, sep=""), size=2.5) +
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
            ylab("TSS Enrich Score") +
            xlab("FRiP") + 
            ggtitle("") 






### MERGE METADATA ATAC-SEQ WITH WGS ---
df <- merge(metadata, metadata_wgs, by.x="Sample_ID", by.y="SampleID")
r4 <- getCorr(x=df$Purity, y=df$ATACseq_FRiP)

#> cor.test(x=df$Purity, y=df$ATACseq_FRiP, method="spearman", na.action = "na.exclude")
#
#        Spearman's rank correlation rho
#
#data:  df$Purity and df$ATACseq_FRiP
#S = 27540, p-value = 4.378e-05
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#      rho
#0.4743584
#
#Warning message:
#In cor.test.default(x = df$Purity, y = df$ATACseq_FRiP, method = "spearman",  :
#  Cannot compute exact p-value with ties


### PLOT ---
p4 <- ggplot(df, aes(x=Purity, y=ATACseq_FRiP)) + 
            stat_smooth(method="lm", geom = "smooth", formula = y ~ x, position = "identity", size=0.5, fullrange = FALSE, se=TRUE, na.rm=TRUE) +	
            geom_point(fill="#e41a1c", color="#e41a1c", shape = 21, stroke=0.5, size=0.8, alpha=0.7, na.rm=TRUE) +
            coord_cartesian(xlim=c(0.2, 1), ylim=c(0, 0.3)) +
            scale_x_continuous(breaks=seq(0.2,1, by=0.2), labels=seq(0.2,1, by=0.2)) +
            scale_y_continuous(breaks=seq(0,0.3, by=0.05), labels=seq(0,0.3, by=0.05) ) +
            annotate("text", x = 0.4, y = 0.25, label = paste("R = ", r4, sep=""), size=2.5) +	
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
            ylab("FRiP") +
            xlab("Purity") + 
            ggtitle("")


### GENERATE PLOT ---
file.plot <- file.path(dir.reproduce_fig, "supplementary_figure_03_b_e.pdf")
pdf(file.plot, width=6.5, height=2)
    grid.arrange(p1, p2, p3, p4, nrow=1, ncol=4)
dev.off()




