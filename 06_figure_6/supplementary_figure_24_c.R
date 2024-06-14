###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### LOAD LIBRARIES ---
library("stringr")
library("data.table")
library("motifStack")
library("TFBSTools")

### DEFINE PATH ---
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/analysis/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

### DEFINE FILES ---
file.rds_motifs <- file.path(dir.reproduce, "reference/jaspar2018_plus_AR_pfm.rds")


### GENES ---
tfs <- c("ZNF263","AR","MYC")

### MOTIF ---
motif_prefix_znf <- "MA0528.1_ZNF263"
motif_prefix_myc <- "MA0147.3_MYC"
motif_prefix_ar <- "ExtendedSite_AR"
motif_prefix <- c(motif_prefix_znf, motif_prefix_ar, motif_prefix_myc)

### PREPARE DESIGN TABLE ---
dat <- data.frame(motif_prefix=motif_prefix, gene=tfs)
dat$MotifID <- unlist(lapply(stringr::str_split(dat$motif_prefix, "_"), function(x) x[1]))

### LOAD MOTIFS ---
list.tfs <- readRDS(file=file.rds_motifs)



### EXTRACT MOTIFS ---
list.motif_tf_pfm <- list()
#list.motif_tf_pcm <- list()
for(i in 1:nrow(dat)){
    gene <- dat$gene[i]
    motif_id <- dat$MotifID[which(dat$gene == gene)]
    
    cat("START:", gene, motif_id, "\n", sep=" ")

    mat_pfm <- motifStack::pcm2pfm(as.data.frame(list.tfs[[motif_id]]@profileMatrix))
    #mat_pcm <- as.data.frame(list.tfs[[motif_id]]@profileMatrix)

    list.motif_tf_pfm[[gene]] <- methods::new("pfm", mat=as.matrix(mat_pfm), name=gene)
    #list.motif_tf_pcm[[gene]] <- methods::new("pcm", mat=as.matrix(mat_pcm), name=gene)

    cat("PROCESSED:", gene, motif_id, "\n", sep=" ")
}


### PLOT MOTIF LOGO STACK ---
file.plot <- file.path(dir.reproduce_fig, "supplementary_figure_24_c.pdf")
pdf(file.plot, width=4, height=4)
    motifStack::motifStack(pfms=list.motif_tf_pfm, layout="stack")
dev.off()


