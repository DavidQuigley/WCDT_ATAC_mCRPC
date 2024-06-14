###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### LOAD LIBRARIES ---
library("stringr")
library("data.table")
library("reshape2")
library("ggplot2")
library("gridExtra")

### DEFINE PATH ---
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/analysis/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

dir.footprints <- file.path(dir.reproduce_data, "tf_footprints/mcrpc_subtypes")


### DEFINE FILES ---
file.rds_footprints <- file.path(dir.footprints, "data_footprints/tobias_mcrpc_subtypes_tf_footprints_featureids_bound_unbound.rds")
file.rds_tf_hits <- file.path(dir.footprints, "data_footprints/tobias_mcrpc_subtypes_tf_hits.rds")

### LOAD FOOTPRINTS ---
list.rds_footprints <- readRDS(file=file.rds_footprints)

### LOAD TF HITS ---
list.rds_tf_hits <- readRDS(file=file.rds_tf_hits)
tfs <- list.rds_tf_hits$tf_footprints_hits_score$motif_prefix

### GET BOUND/UNBOUND FREQUENCY AND PERCENTAGE ---
list.dat <- list()
for(motif_prefix in tfs){

    freq_bound <- length( unique( unlist(list.rds_footprints$bound[[motif_prefix]]) ) )
    freq_unbound <- length( unique( unlist(list.rds_footprints$unbound[[motif_prefix]]) ) )
    total_footprints <- freq_bound + freq_unbound
    percent_bound <- round( ( freq_bound / total_footprints )*100, 3)
    percent_unbound <- round( ( freq_unbound / total_footprints )*100, 3)

    list.dat[[motif_prefix]] <- data.frame(motif=motif_prefix, freq_bound, freq_unbound, total_footprints, percent_bound, percent_unbound)
    
}

### AGGREGATE DATA ---
dat <- do.call(rbind.data.frame, list.dat)
rownames(dat) <- NULL



##########################################################################################################
### GET DATA FOR PERCENT ---
dat_frac <- subset(dat, select=c("motif","percent_bound","percent_unbound"))
dat_frac <- dat_frac[order(dat_frac$percent_bound, decreasing=FALSE),]

### SEPARATE DATA ---
dat_frac1 <- dat_frac[102:203,]
dat_frac2 <- dat_frac[1:101,]

dm_frac1 <- reshape2::melt(dat_frac1, id.vars="motif", variable.name="item", value.name="percent")
dm_frac2 <- reshape2::melt(dat_frac2, id.vars="motif", variable.name="item", value.name="percent")

### FACTORIZE ---
dm_frac1$motif <- factor(dm_frac1$motif, levels=dat_frac1$motif)
dm_frac2$motif <- factor(dm_frac2$motif, levels=dat_frac2$motif)

dm_frac1$item <- factor(dm_frac1$item, levels=c("percent_unbound","percent_bound"))
dm_frac2$item <- factor(dm_frac2$item, levels=c("percent_unbound","percent_bound"))



### FUNCTION: getStackedBarplot() ---
getStackedBarplot <- function(dm){
    # PLOT ---
    p <- ggplot(dm, aes(x=percent, y=motif)) +
            geom_bar(aes(fill=item), stat="identity", color="#000000", width=0.8, size=0.2) +
            scale_fill_manual(values=c("#969696","#000000")) +
            coord_cartesian(xlim=c(0,100)) +
            scale_x_continuous(breaks=seq(0,100, by=20), labels=seq(0,100, by=20)) +
		    theme(
			    axis.text.x = element_text(size = 4.5, color="#000000"),
			    axis.text.y = element_text(size = 4.5, color="#000000"),
			    axis.title = element_text(size = 5, color="#000000"),
			    plot.title = element_text(size = 5, color="#000000", hjust=0.5),
			    panel.grid.major.y = element_blank(),
                panel.grid.major.x = element_line(size=0.1, color="#BDBDBD"),
			    panel.grid.minor = element_blank(),
			    axis.ticks = element_line(size = 0.2, color="#000000"),	
			    panel.background = element_rect(fill = "#FFFFFF", color = "#000000"),
			    legend.text = element_text(size = 5, color="#000000"),
			    legend.title = element_blank(),
			    legend.key.size = unit(0.2, "cm"),			
			    legend.position="none",
                legend.box = "horizontal") + 
            guides(fill = guide_legend(nrow = 2)) +                   
		    ylab("") +            
		    xlab("% of TF footprints") + 
            ggtitle("") 

    return(p)
}



### GET PLOT ---
p1 <- getStackedBarplot(dm=dm_frac1)
p2 <- getStackedBarplot(dm=dm_frac2)

### GENERATE PLOT ---
file.plot <- file.path(dir.reproduce_fig, "supplementary_figure_14.pdf")
pdf(file.plot, width=6.8, height=8)
    grid.arrange(p1, p2, nrow=1, ncol=2)
dev.off()


