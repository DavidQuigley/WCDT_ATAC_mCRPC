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
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/analysis/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

dir.footprints <- file.path(dir.reproduce_data, "tf_footprints/mcrpc_subtypes")

### DEFINE FILES ---
file.rds_bindetect_results <- file.path(dir.footprints, "data_footprints/tobias_mcrpc_subtypes_bindetect_results.rds")

### MOTIF ---
motif_prefix <- "MA0528.1_ZNF263"

### GROUPS ---
subtypes <- c("ARpNEn","ARpNEp","ARlNEn","ARnNEp","ARnNEn")



### GET PAIRWISE ALL GROUP COMBINATIONS ---
dm_grp <- do.call(rbind.data.frame, combn(x=subtypes, m=2, simplify=FALSE)) %>%
                dplyr::mutate_all(as.character)
colnames(dm_grp) <- c("Group1","Group2")
rownames(dm_grp) <- NULL
dm_grp$id <- apply(dm_grp, 1, function(x) paste(x[1], x[2], sep="_"))
dm_grp$f_name <- apply(dm_grp, 1, function(x) paste("tobias_diff", paste(x[1], x[2], sep="_"), sep="_"))



### FUNCTION: parseData() ---
parseData <- function(dat, id){
    # GET GROUPS ---
    grp1 <- stringr::str_split(id, "_")[[1]][1]
    grp2 <- stringr::str_split(id, "_")[[1]][2]

    # CHANGE COLNAMES ---
    colnames(dat)[stringr::str_detect(colnames(dat), "change")] <- "change"
    colnames(dat)[stringr::str_detect(colnames(dat), "pvalue")] <- "pvalue"

    colnames(dat)[which(colnames(dat) == paste(grp1, "mean_score", sep="_"))] <- "group1_mean_score"
    colnames(dat)[which(colnames(dat) == paste(grp2, "mean_score", sep="_"))] <- "group2_mean_score"

    colnames(dat)[which(colnames(dat) == paste(grp1, "bound", sep="_"))] <- "group1_bound"
    colnames(dat)[which(colnames(dat) == paste(grp2, "bound", sep="_"))] <- "group2_bound"

    # ADD ID ---
    dat$ID <- id

    return(dat)
}



######################################################################
### LOAD TOBIAS BINDETECT RESULTS ---
list.rds_bindetect_results <- readRDS(file=file.rds_bindetect_results)


### LOAD DATA BY COMPARISON ID ----
list.dat <- list()
for(i in 1:nrow(dm_grp)){
    id <- dm_grp$id[i]

    cat("START:", id, "\n", sep=" ")

    # LOAD DATA ---
    dat <- list.rds_bindetect_results[[id]]

    # PARSE DATA ---
    list.dat[[id]] <- parseData(dat, id)

    cat("DONE:", id, "\n", sep=" ")
}



### EXTRACT DATA BY MOTIF ---
list.dat_znf <- lapply(list.dat, function(x) subset(x, x$output_prefix == motif_prefix))

### AGGREGATE DATA ---
dat <- do.call(rbind.data.frame, list.dat_znf)
rownames(dat) <- NULL

### TRIM DATA ---
df <- subset(dat, select=c("ID","group1_bound","group2_bound"))
colnames(df) <- c("ID","group1","group2")

### MELT DATA ---
dm <- reshape2::melt(df, id.vars="ID", variable.name="group", value.name="value")

### FACTORIZE ---
dm$ID <- factor(dm$ID, levels=dat$ID)
dm$group <- factor(dm$group, levels=c("group1","group2"))


### PLOT ---
cpallete <- c("#FFFFFF","#FFFFFF")
p <- ggplot(dm, aes(x=group, y=value, fill=group)) + 
		geom_bar(stat="identity", color="#000000", width=0.8, size=0.2) + 
		scale_fill_manual(values=cpallete) +
        facet_wrap(~ID, nrow=2, ncol=5, scales="free", strip.position="top") +
		theme(
			axis.text.x = element_text(size = 6, color="#000000"),
			axis.text.y = element_text(size = 6, color="#000000"),
			axis.title = element_text(size = 8, color="#000000"),
			plot.title = element_text(size = 10, color="#000000", hjust=0),
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			axis.ticks = element_line(size = 0.2, colour="#000000"),	
            strip.text = element_text(size=5, color="#000000"),
            strip.background = element_rect(fill="#FFFFFF", color="#FFFFFF"),
			panel.background = element_rect(fill = "#FFFFFF", colour = "#000000"),           
			legend.text = element_text(size = 5, color="#000000"),
			legend.title = element_blank(),
			legend.key.size = unit(0.3, "cm"),			
			legend.position="none") + 
			ylab("") +
			xlab("") + 
            ggtitle("") 

### GENERATE PLOT ---
file.plot <- file.path(dir.reproduce_fig, "supplementary_figure_22_e.pdf")
pdf(file.plot, height=3, width=6)
    grid.arrange(p, nrow=1, ncol=1)
dev.off()


