###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### LOAD LIBRARIES ---
library("stringr")
library("data.table")
library("RColorBrewer")
library("circlize")
library("ComplexHeatmap")

### DEFINE PATH ---
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/analysis/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

dir.footprints <- file.path(dir.reproduce_data, "tf_footprints/mcrpc_subtypes")

### DEFINE FILES ---
file.rds_tf_hits <- file.path(dir.footprints, "data_footprints/tobias_mcrpc_subtypes_tf_hits.rds")

### LOAD TF HITS ---
list.rds_tf_hits <- readRDS(file=file.rds_tf_hits)

### GET META SCORE DATA ---
mat <- list.rds_tf_hits$tf_footprints_hits_score
rownames(mat) <- mat$gene
mat$motif_prefix <- NULL
mat$gene <- NULL



### NORMALIZE DATA ---
mat_norm <- mat
for(i in 1:nrow(mat)){
	mat_norm[i,] <- mat[i,] / min(setdiff(mat[i,],0))
}


### GENES TO HIGHLIGHT -----
genes_highlight <- c("AR","FOXA1","HOXB13","ASCL1","NEUROD1","SP1","SP2","KLF5","POU3F2","ZNF263")

index_genes_highlight <- which(rownames(mat) %in% genes_highlight)

text_size <- rep(0.29, nrow(mat))
text_size[index_genes_highlight] <- 0.4

text_color <- rep("#000000", nrow(mat))
text_color[index_genes_highlight] <- "#e31a1c"


### GENERATE COLOR PALETTE ---
jColFun <- colorRampPalette(brewer.pal(n = 9, "Reds"))
col_fun1 <- circlize::colorRamp2(breaks=c(0:5), colors=jColFun(6), transparency=0, space="RGB")

### PLOT ---
file.plot <- file.path(dir.reproduce_fig, "figure_4_e.pdf")
pdf(file.plot, width=4, height=4)
    circlize::circos.clear()
    circlize::circos.par(gap.after = c(10))
    circlize::circos.heatmap(mat=mat_norm, split=NULL, col=col_fun1, na.col = "#DDDDDD", 
                            cell.border="#000000", cell.lty=1, cell.lwd=0.1, 
                            bg.border=NA, bg.lty=NA, bg.lwd=NA, ignore.white=TRUE, 
                            cluster=TRUE, clustering.method="ward.D2", distance.method="euclidean", 
                            rownames.side="outside", rownames.cex=text_size, rownames.col=text_color, 
                            show.sector.labels=FALSE, cell_width=5)

    circlize::circos.clear()

    lgd <- ComplexHeatmap::Legend(col_fun = col_fun1)
    grid.draw(lgd)

dev.off()


