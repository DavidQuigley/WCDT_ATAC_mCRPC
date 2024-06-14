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

### DEFINE FILES ---
file.dat <- file.path(dir.reproduce_data, "varbin_correlation_atacseq_rnaseq.tsv.gz")

### ATAC-SEQ - RNA-SEQ CORR ---
mat <- data.table::fread(file.dat, data.table=FALSE, header=TRUE, stringsAsFactors=FALSE)
rownames(mat) <- mat$V1
mat$V1 <- NULL
mat <- na.omit(mat)

### GET MAX. CORRELATION ---
dat <- tibble::tibble(data.frame(EnsemblID=rownames(mat), Corr=apply(mat, 1, function(x) max(x))))

### GET MEAN CORRELATION ---
y_mean <- round(mean(dat$Corr), 2)



### PLOT ---
p <- ggplot(dat, aes(x=Corr)) +
        geom_vline(xintercept=y_mean, color="gray70", size=0.4, linetype=2) +
        geom_density(fill="#e41a1c", color="#000000", alpha=0.5, na.rm=TRUE) +
        scale_y_continuous(breaks=seq(0, 3, by=0.5), labels=seq(0, 3, by=0.5)) +
        coord_cartesian(xlim=c(-1,1), ylim=c(0,3)) +
        theme(
            axis.text.x = element_text(size = 6, color="#000000"),
            axis.text.y = element_text(size = 6, color="#000000"),
            axis.title = element_text(size = 6, color="#000000"),
            plot.title = element_text(size = 6, color="#000000", hjust=0.5),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.ticks = element_line(size=0.2, color="#000000"),    
            strip.text = element_text(size=6, color="#000000"),
            strip.background = element_rect(fill="#FFFFFF", color="#FFFFFF"),
            panel.background = element_rect(fill="#FFFFFF", color="#000000"),
            legend.text = element_text(size = 6, color="#000000"),
            legend.title = element_blank(),
            legend.key.size = unit(0.2, "cm"),
            legend.position = "none") +
        ylab("Density") +
        xlab("Spearman's Correlation Coefficient (R)") + 
        ggtitle("") 


### GENERATE PLOT ---
file.plot <- file.path(dir.reproduce_fig, "supplementary_figure_07_b.pdf")
pdf(file.plot, width=2, height=2)
    grid.arrange(p, nrow=1, ncol=1)
dev.off()


