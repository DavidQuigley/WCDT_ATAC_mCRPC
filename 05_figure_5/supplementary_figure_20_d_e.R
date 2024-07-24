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

### GET PEAK AND GENE DEGREE ---
df_peak_degree <- dat %>% dplyr::count(Peak, sort=TRUE, name="Freq")
df_gene_degree <- dat %>% dplyr::count(Gene, sort=TRUE, name="Freq")


### PLOT: PEAK DEGREE HISTOGRAM ---
p1 <- ggplot(df_peak_degree, aes(x=Freq)) +
            geom_histogram(stat="count", fill="#000000", color="#000000") +
            coord_cartesian(xlim=c(0,30), ylim=c(0,15000) ) +
            scale_x_continuous(breaks=seq(0,30, by=5) ) +
            scale_y_continuous(breaks=seq(0,15000, by=2500) ) +
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
            ylab("No. of ATAC-seq peaks") +
            xlab("No. of genes mapped per ATAC-seq peaks") + 
            ggtitle("")


### PLOT: GENE DEGREE HISTOGRAM ---
p2 <- ggplot(df_gene_degree, aes(x=Freq)) +
            geom_histogram(stat="count", fill="#000000", color="#000000") +
            coord_cartesian(xlim=c(0,120), ylim=c(0,3000) ) +
            scale_x_continuous(breaks=seq(0,120, by=30) ) +
            scale_y_continuous(breaks=seq(0,3000, by=500) ) +
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
            ylab("No. of genes") +
            xlab("No. of ATAC-seq peaks linked per gene") + 
            ggtitle("")

### GENERATE PLOT ---
file.plot <- file.path(dir.reproduce_fig, "supplementary_figure_20_d_e.pdf")
pdf(file.plot, height=2, width=4)
    grid.arrange(p1, p2, nrow=1, ncol=2)
dev.off()

  

