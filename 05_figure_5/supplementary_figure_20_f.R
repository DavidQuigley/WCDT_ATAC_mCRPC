###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### LOAD LIBRARIES ---
library("stringr")
library("data.table")
library("ggplot2")
library("gridExtra")

### DEFINE PATH ---
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

dir.tfpeaks2gene <- file.path(dir.reproduce_data, "tfpeaks2gene")

### DEFINE FILES ---
file.dat <- file.path(dir.tfpeaks2gene, "wcdt_tfpeaks2genes_correlated_r_0p4_pval_0p05_distTSS_500kb.tsv.gz")

### LOAD PEAKS2GENE DATA ---
dat <- data.table::fread(file=file.dat, sep="\t", header=TRUE, nThread=60, data.table=FALSE, verbose=FALSE)

### GET POSITION GROUP ---
dat$distTSS <- round(dat$distTSS/1000, 2)
dat$Group <- cut(dat$distTSS, seq(-500,500, by=10), include.lowest = TRUE)


### FUNCTION: plotCorrDistribution() ---
plotCorrDistribution <- function(dat){
    p <- ggplot(dat, aes(x=Group)) +
            geom_histogram(stat="count", fill="#006094", color=NA) +
            coord_cartesian(ylim=c(0,2500) ) +
            scale_y_continuous(breaks=seq(0,2500, by=500) ) +
            scale_x_discrete(breaks=c("[-500,-490]", "(0,10]", "(490,500]"), labels=c("-500","TSS","500")) +
            theme(
               axis.text.x = element_text(size = 5, color="#000000"),
               axis.text.y = element_text(size = 5, color="#000000"),
               axis.title = element_text(size = 5, color="#000000"),
               plot.title = element_text(size = 5, color="#000000", hjust=0.5),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),  
               axis.ticks.x = element_line(size=0.2, color="#000000"),  
               axis.ticks.y = element_line(size=0.2, color="#000000"),
               strip.text = element_text(size=5, color="#000000"),
               strip.background = element_rect(fill="#FFFFFF", color="#FFFFFF"),
               panel.background = element_rect(fill="#FFFFFF", color="#000000"),
               legend.text = element_text(size = 5, color="#000000"),
               legend.title = element_blank(),
               legend.key.size = unit(0.3, "cm"),
               legend.position = "none") +
            ylab("No. of highly correlated \n peak-gene pairs") +
            xlab("Distance to TSS (kbp)") + 
            ggtitle("")
    
    return(p) 
}

# GET PLOT ---
p <- plotCorrDistribution(dat)

### GENERATE PLOT ---
file.plot <- file.path(dir.reproduce_fig, "supplementary_figure_20_f.pdf")
pdf(file.plot, width=2, height=2)
        grid.arrange(p, nrow=1, ncol=1)
dev.off()
