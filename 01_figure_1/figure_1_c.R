###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### DEFINE LIBRARIES ---
library("stringr")
library("data.table")

### DEFINE PATH ---
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

dir.bigwig <- file.path(dir.wrk, "secondary_data/wcdt_atacseq/bigwig_merge")

### DEFINE FILES ----
file.rds_atacseq_chromatin_variants <- file.path(dir.reproduce_data, "atacseq_chromatin_variants_benign_localizedpca_mcrpc_adeno_nepc.rds")
file.utility_functions <- file.path(dir.reproduce_scripts, "00_prepare_data/00_01_utility_functions.R")

### LOAD FUNCTIONS ---
source(file.utility_functions)

### CHROMOSOMES ---
chromosomes <- paste("chr", c(1:22,"X","Y"), sep="")

### FEATURE ---
stages <- c("benign","pca","mcrpc_adeno","mcrpc_nepc")
stages_bw <- c("pome_benign","tcga_pca","mcrpc_adeno","ARnNEp")

### LOAD ATAC-SEQ CHROMATIN VARIANTS DATA ---
list.rds_atacseq_chromatin_variants <- readRDS(file=file.rds_atacseq_chromatin_variants)

### GET CHROMATIN VARIANT DATA BY FEATURES-STAGE ---
dat <- list.rds_atacseq_chromatin_variants$atacseq_diffpeaks_up_open
dat <- dat[which(dat$GeneType == "protein_coding"),]


### FUNCTION: getBED() ---
getBED <- function(feature_id){
    list.feat <- stringr::str_split(feature_id, "_")
    chr <- unlist(lapply(list.feat, function(x) x[1]))
    start <- as.numeric(unlist(lapply(list.feat, function(x) x[2])))
    end <- as.numeric(unlist(lapply(list.feat, function(x) x[3])))

    # GET BED DATA ---
    bed <- data.frame(chr=chr, start=start, end=end)

    # GET GRANGES ---
    gr_bed <- GenomicRanges::makeGRangesFromDataFrame(bed, keep.extra.columns=FALSE)

    return(gr_bed)
}

### EXTRACT BED FILES ---
list.bed <- list()
for(stage in stages){
    index_stage <- which(colnames(dat) == stage)
    feature_id <- dat$FeatureID[which(dat[,index_stage] == 1)]
    
    # GET BED ---
    list.bed[[stage]] <- getBED(feature_id)
}




### TILE GENOME ---
gr_bins <- binGenome(genome=BSgenome.Hsapiens.UCSC.hg38, chromosome=chromosomes, tile.size=500)

### BINS TO FEATURE_IDS ---
df.bins <- as.data.frame(gr_bins, stringsAsFactors=FALSE)
dt.bins <- data.table::as.data.table(df.bins)

### GET FEATURE IDS ---
mod_cols <- c("seqnames","start","end") 
myfun <- function(y) paste(y$seqnames, y$start, y$end, sep = "_")
feature_ids <- dt.bins[, myfun(.SD), .SDcols = mod_cols]
gr_bins$feature_ids <- feature_ids


### GET GENOME BIN TO BEN DATA OVERLAP ----
list.gr_bed_bins <- list()
for(stage in stages){
    # GET OVERLAPS ---
    list.gr_bed_bins[[stage]] <- IRanges::subsetByOverlaps(x=gr_bins, ranges=list.bed[[stage]])
}


### GET FEATURE IDS ---
fids_1 <- list.gr_bed_bins[["benign"]]$feature_ids
fids_2 <- list.gr_bed_bins[["pca"]]$feature_ids
fids_3 <- list.gr_bed_bins[["mcrpc_adeno"]]$feature_ids
fids_4 <- list.gr_bed_bins[["mcrpc_nepc"]]$feature_ids

### GET PEAKS BIN IDENTIFIED IN EXCLUSIVELY IN THE STAGE ---
fids_only_1 <- setdiff(fids_1, unique(c(fids_2, fids_3, fids_4)) )
fids_only_2 <- setdiff(fids_2, unique(c(fids_1, fids_3, fids_4)) )
fids_only_3 <- setdiff(fids_3, unique(c(fids_2, fids_1, fids_4)) )
fids_only_4 <- setdiff(fids_4, unique(c(fids_2, fids_3, fids_1)) )


### REFINED BED REGION -----
list.gr_bed_bins_refined <- list()
list.gr_bed_bins_refined[["benign"]] <- subset(list.gr_bed_bins[["benign"]], list.gr_bed_bins[["benign"]]$feature_ids %in% fids_only_1)
list.gr_bed_bins_refined[["pca"]] <- subset(list.gr_bed_bins[["pca"]], list.gr_bed_bins[["pca"]]$feature_ids %in% fids_only_2)
list.gr_bed_bins_refined[["mcrpc_adeno"]] <- subset(list.gr_bed_bins[["mcrpc_adeno"]], list.gr_bed_bins[["mcrpc_adeno"]]$feature_ids %in% fids_only_3)
list.gr_bed_bins_refined[["mcrpc_nepc"]] <- subset(list.gr_bed_bins[["mcrpc_nepc"]], list.gr_bed_bins[["mcrpc_nepc"]]$feature_ids %in% fids_only_4)


### WRITE BED FILES ---
list.files_bed <- list()
list.files_matrix <- list()
for(stage in stages){
    # GET BED ---
    bed <- list.gr_bed_bins_refined[[stage]]

    # WRITE OUTPUT BED FILE ---
    list.files_matrix[[stage]] <- tempfile(pattern=paste("diffpeaks", stage, sep="_"), tmpdir=tempdir(), fileext=".tsv.gz")
    list.files_bed[[stage]] <- tempfile(pattern=paste("diffpeaks", stage, sep="_"), tmpdir=tempdir(), fileext=".bed")
    write.table(bed, list.files_bed[[stage]], sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
}


### GET FILES ---
files.bigwig <- sapply(stages_bw, function(x) { paste(dir.bigwig, "/", sprintf("%s", x), ".bw", sep="") })
files.bigwig <- paste(files.bigwig, collapse=" ")



### FUNCTION: getComuteMatrix() ---
getComuteMatrix <- function(files.bigwig, files.bed, file.matrix, sampleids, n_cores=10){
    # 'computeMatrix' FUNCTION FROM 'deepTools'
    # https://deeptools.readthedocs.io/en/latest/index.html

    cmd <- paste("computeMatrix scale-regions",
                "--scoreFileName", files.bigwig,
                "--regionsFileName", files.bed,
                "--outFileName", file.matrix,
                "--regionBodyLength", 500,
                "--beforeRegionStartLength", 1000,
                "--afterRegionStartLength", 1000,
                "--binSize", 5,
                "--missingDataAsZero",
                "--samplesLabel", sampleids,
                "--verbose",
                "-p", n_cores,  
                sep=" ")
    
    return(cmd)
}


### FUNCTION: getPlotHeatmap() ---
getPlotHeatmap <- function(file.matrix, file.plot, sample_index){
    # 'plotHeatmap' FUNCTION FROM 'deepTools'
    # https://deeptools.readthedocs.io/en/latest/index.html

    cmd <- paste("plotHeatmap",
                "--matrixFile", file.matrix,
                "--outFileName", file.plot,
                "--sortRegions", "descend",
                "--sortUsing", "mean",
                "--sortUsingSamples", sample_index,
                "--colorMap", "Reds",
                "--interpolationMethod", "gaussian",
                #"--whatToShow", "\"heatmap only\"",
                "--whatToShow", "\"heatmap and colorbar\"",
                #"--whatToShow", "\"plot and heatmap\"",
                #"--whatToShow", "\"plot, heatmap and colorbar\"",
                "--startLabel", "\"Start\"",
                "--endLabel", "\"End\"",
                "--regionsLabel bed_peaks", 
                "--samplesLabel benign pca mcrpc nepc",
                "--legendLocation", "\"best\"",
                #"--regionsLabel", "\"ATAC-seq Peaks\"",
                "--xAxisLabel", "\"Peaks\"",
                "--linesAtTickMarks", 
                "--heatmapHeight 2",
                "--heatmapWidth 2",
                #"--boxAroundHeatmaps no",            
                #"--kmeans 2",                
                "--dpi", 600,
                "--verbose",
                sep=" ")
    
    return(cmd)
}


### GET COMPUTE MATRIX ---
for(stage in stages){
    cmd <- getComuteMatrix(files.bigwig=files.bigwig, 
                            files.bed=list.files_bed[[stage]], 
                            file.matrix=list.files_matrix[[stage]], 
                            sampleids=paste(stages, collapse=" "), 
                            n_cores=60)
    system(cmd)
}


### SAMPLE INDICES ---
list.sample_index <- c("1 2 3 4", "2 1 3 4", "3 4 2 1", "4 3 2 1")
names(list.sample_index) <- stages

### PLOT HEATMAP ---
list.files_plot <- list()
for(stage in stages){
    list.files_plot[[stage]] <- file.path(dir.reproduce_fig, paste("figure_1_c_", stage, ".pdf", sep=""))

    cmd <- getPlotHeatmap(file.matrix=list.files_matrix[[stage]], 
                            file.plot=list.files_plot[[stage]], 
                            sample_index=list.sample_index[[stage]])
    system(cmd)
}

### MERGE PLOT PDFs ---------------------------------------------------
file.merge_plot <- file.path(dir.reproduce_fig, "figure_1_c.pdf")
cmd_qpdf <- paste("qpdf --empty --pages",
                paste(unlist(list.files_plot), collapse=" "),
                "--", file.merge_plot,
                sep=" ")

system(cmd_qpdf)

cat("\n", "FILE GENERATED:", file.merge_plot, "\n", sep=" ")
