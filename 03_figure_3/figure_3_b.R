###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### LOAD LIBRARIES ---
library("stringr")
library("data.table")
library("GenomicRanges")
library("IRanges")

### DEFINE PATH ---
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/analysis/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

dir.bedgraph <- file.path(dir.wrk, "secondary_data/wcdt_atacseq/bedgraph_merge")
dir.peaks <- file.path(dir.wrk, "secondary_data/wcdt_atacseq/peaks/subtypesAR")

### DEFINE FILES ----
file.rds_atacseq_masterdata <- file.path(dir.reproduce_data, "wcdt_atacseq_masterdata.rds")
file.rds_chipseq_bed <- file.path(dir.reproduce_data, "chipseq_bed.rds")
file.kwtest <- file.path(dir.reproduce_data, "wcdt_atacseq_mcrpc_subtypes_KruskalWallisTest_summary.tsv.gz")
file.h3k27ac_hichip <- file.path(dir.reproduce_data, "hg38_h3k27ac_lncap_hichip_giambartolomei_etal_ajhg_2021_interactions.bedpe")
file.utility_functions <- file.path(dir.reproduce_scripts, "00_prepare_data/00_01_utility_functions.R")
file.plot_functions <- file.path(dir.reproduce_scripts, "00_prepare_data/00_02_plot_tracks_functions.R")

file.gene_annotation <- file.path(dir.wrk, "reference/annotation_data/gene_annotation/gencode_v28/gencode.v28.annotation.gtf.gz")

### LOAD FUNCTIONS ---
source(file.utility_functions)
source(file.plot_functions)

### SUBTYPES ---
subtypes <- c("ARpNEn","ARlNEn","ARpNEp","ARnNEp","ARnNEn")

### COLOR ---
cp_subtypes <- c("#e31a1c","#fb9a99","#ff7f00","#33a02c","#1f78b4")



### SUBTYPES BEDGRAPH AND ATAC-SEQ PEAK FILES ---
files.bedgraph <- sapply(subtypes, function(x) { paste(dir.bedgraph, "/", sprintf("%s", x), ".bedgraph.gz", sep="") })
files.peaks <- sapply(subtypes, function(x) { paste(dir.peaks, "/subtypesAR_narrow_peaks_", sprintf("%s", x), ".bed", sep="") })



#########################################################################################################

### LOAD CHIP-SEQ DATA ---
list.rds_chipseq_bed <- readRDS(file=file.rds_chipseq_bed)
gr_h3k27ac <- list.rds_chipseq_bed$pomerantz_mcrpc$H3K27ac

### LOAD H3K27AC HICHIP DATA ---
#dat_h3k27ac_hichip <- data.table::fread(file=file.h3k27ac_hichip, sep="\t", header=TRUE, nThread=50, data.table=FALSE, verbose=FALSE)

### GET DIFFPEAKS BED DATA ---
dat_kwtest <- data.table::fread(file=file.kwtest, sep="\t", header=TRUE, nThread=50, data.table=FALSE, verbose=FALSE)
colnames(dat_kwtest)[1] <- "FeatureID"
dat_kwtest <- dat_kwtest[order(dat_kwtest$pvalue, decreasing=FALSE),]
dat_kwtest <- subset(dat_kwtest, dat_kwtest$pvalue <= 0.001)

dat_kwtest$chr <- unlist(lapply(stringr::str_split(dat_kwtest$FeatureID, "_"), function(x) x[1]))
dat_kwtest$start <- as.numeric(unlist(lapply(stringr::str_split(dat_kwtest$FeatureID, "_"), function(x) x[2])))
dat_kwtest$end <- as.numeric(unlist(lapply(stringr::str_split(dat_kwtest$FeatureID, "_"), function(x) x[3])))
diffpeaks_bed <- subset(dat_kwtest, select=c("chr","start","end"))


#########################################################################################################
### GET BEDGRAPH ---
list.bedgraph <- loadBedGraph(files.bedgraph, n_cores=50)

### GET PEAKS ---
list.peaks <- loadPeaks(files.peaks, n_cores=10)

### GET REF ANNOTATION DATA FOR GENE AND EXON ---
annot.gene <- parseGTF(file.gtf=file.gene_annotation, feature.type="gene")
annot.exon <- parseGTF(file.gtf=file.gene_annotation, feature.type="exon")




#########################################################################################################
### PLOT BY GENE: 1 ---
gene <- "AR"
file.plot_1 <- file.path(dir.reproduce_fig, paste("figure_3_b_", gene, ".pdf", sep=""))
getPlotbyGene(annot.gene, annot.exon, list.peaks, list.bedgraph, 
                diffpeaks_bed, gr_h3k27ac, file.plot=file.plot_1, subtypes, cp_subtypes, 
                upstream=800000, downstream=100000, range=c(0,3), height=5, width=2, textyPos=2.5, gene=gene)
     

### PLOT BY GENE: 2 ---
gene <- "NKX3-1"
file.plot_2 <- file.path(dir.reproduce_fig, paste("figure_3_b_", gene, ".pdf", sep=""))
getPlotbyGene(annot.gene, annot.exon, list.peaks, list.bedgraph, 
                diffpeaks_bed, gr_h3k27ac, file.plot=file.plot_2, subtypes, cp_subtypes, 
                upstream=12000, downstream=5000, range=c(0,1.2), height=5, width=2, textyPos=0.7, gene=gene)
     


### PLOT BY GENE: 3 ---
gene <- "GPR37L1"
file.plot_3 <- file.path(dir.reproduce_fig, paste("figure_3_b_", gene, ".pdf", sep=""))
getPlotbyGene(annot.gene, annot.exon, list.peaks, list.bedgraph, 
                diffpeaks_bed, gr_h3k27ac, file.plot=file.plot_3, subtypes, cp_subtypes, 
                upstream=60000, downstream=10000, range=c(0,3), height=4.5, width=2, textyPos=1, gene=gene)




### MERGE PLOT PDFs ---------------------------------------------------
file.merge_plot <- file.path(dir.reproduce_fig, "figure_3_b.pdf")
cmd_qpdf <- paste("qpdf --empty --pages",
                paste(file.plot_1, file.plot_2, file.plot_3, sep=" "),
                "--", file.merge_plot,
                sep=" ")

system(cmd_qpdf)

cat("\n", "FILE GENERATED:", file.merge_plot, "\n", sep=" ")


