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


### PEAKS-GENE CORRELATION DISTRIBUTION ----
p <- ggplot(dat, aes(x=Corr)) +
            geom_histogram(aes(y=..density..), position="identity", fill="#ff7f00", color=NA, bins=100, alpha=0.7) +
            geom_density() +
            #geom_vline(xintercept=mean(dat$Corr), linetype="dashed", size=0.5, colour="#e31a1c") +
            coord_cartesian(xlim=c(-1,1), ylim=c(0,12) ) +
            scale_x_continuous(breaks=seq(-1,1, by=0.2) ) +
            scale_y_continuous(breaks=seq(0,12, by=2) ) +
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
               legend.key.size = unit(0.5, "cm"),
               legend.position = "none") +
            ylab("Density") +
            xlab("Spearman's Correlation Coefficient (R)") + 
            ggtitle("")
    

### GENERATE PLOT ---
file.plot <- file.path(dir.reproduce_fig, "supplementary_figure_20_b.pdf")
pdf(file.plot, width=2, height=2)
        grid.arrange(p, nrow=1, ncol=1)
dev.off()

