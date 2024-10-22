###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### LOAD LIBRARIES ---
library("stringr")
library("data.table")

### DEFINE PATH ---
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

dir.bigwig <- file.path(dir.wrk, "secondary_data/wcdt_atacseq/bigwig_merge") # REPLACE PATH
dir.data <- tempdir() # REPLACE PATH
dir.scripts <- tempdir() # REPLACE PATH
dir.logs <- tempdir() # REPLACE PATH

### DEFINE FILES ----
file.rds_atacseq_masterdata <- file.path(dir.reproduce_data, "wcdt_atacseq_masterdata.rds")
file.rds_atacseq_read_counts_mcrpc <- file.path(dir.reproduce_data, "atacseq_read_counts_mcrpc.rds")
file.kwtest <- file.path(dir.reproduce_data, "wcdt_atacseq_mcrpc_subtypes_KruskalWallisTest_summary.tsv.gz")

### FEATURES ---
features <- c("promoter","intron","intergenic")
subtypes <- c("ARpNEn","ARlNEn","ARpNEp","ARnNEp","ARnNEn")



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


### FUNCTION: getPlotHeatmapTotal() ---
getPlotHeatmap <- function(file.matrix, file.plot){
    # 'plotHeatmap' FUNCTION FROM 'deepTools'
    # https://deeptools.readthedocs.io/en/latest/index.html

    cmd <- paste("plotHeatmap",
                "--matrixFile", file.matrix,
                "--outFileName", file.plot,
                "--sortRegions", "descend",
                "--sortUsing", "mean",
                "--sortUsingSamples", "5 1 2",
                "--colorMap", "Reds",
                "--interpolationMethod", "gaussian",
                "--whatToShow", "\"heatmap only\"",
                #"--whatToShow", "\"heatmap and colorbar\"",
                #"--whatToShow", "\"plot and heatmap\"",
                #"--whatToShow", "\"plot, heatmap and colorbar\"",
                "--startLabel", "\"Start\"",
                "--endLabel", "\"End\"",
                "--regionsLabel promoter intron intergenic", 
                "--samplesLabel ARpNEn ARlNEn ARpNEp ARnNEp ARnNEn",
                "--legendLocation", "\"best\"",
                #"--regionsLabel", "\"ATAC-seq Peaks\"",
                "--xAxisLabel", "\"Peaks\"",
                "--linesAtTickMarks", 
                "--heatmapHeight 5",
                "--heatmapWidth 2",
                #"--zMin 0.05",
                "--zMax 0.3 0.3 0.3 0.3 0.19",
                #"--boxAroundHeatmaps no",            
                #"--kmeans 2",                
                "--dpi", 600,
                "--verbose",
                sep=" ")
    
    return(cmd)
}


#############################################################################################################################
### LOAD ATAC-SEQ PEAK ANNOTATION DATA ---
list.rds_atacseq_read_counts_mcrpc <- readRDS(file=file.rds_atacseq_read_counts_mcrpc)
annot <- list.rds_atacseq_read_counts_mcrpc$feature_annotation

### LOAD KRUSKAL-WALLIS TEST DATA ---
dat_kwtest <- data.table::fread(file=file.kwtest, sep="\t", header=TRUE, nThread=50, data.table=FALSE, verbose=FALSE)
colnames(dat_kwtest)[1] <- "FeatureID"
dat_kwtest <- merge(dat_kwtest, annot, by="FeatureID")
dat_kwtest <- dat_kwtest[order(dat_kwtest$pvalue, decreasing=FALSE),]
dat_kwtest <- subset(dat_kwtest, dat_kwtest$pvalue <= 0.001)

dat_kwtest <- subset(dat_kwtest, dat_kwtest$Feature %in% c("Promoter","Intron","Distal Intergenic"))
dat_kwtest$Feature[which(dat_kwtest$Feature == "Promoter")] <- "promoter"
dat_kwtest$Feature[which(dat_kwtest$Feature == "Intron")] <- "intron"
dat_kwtest$Feature[which(dat_kwtest$Feature == "Distal Intergenic")] <- "intergenic"
rownames(dat_kwtest) <- NULL


### EXTRACT BED FILES ---
list.beds <- list()
for(feat in features){
    # GET BED ---    
    d <- dat_kwtest[which(dat_kwtest$Feature == feat),]
    bed <- subset(d, select=c("seqnames","start","end"))

    # WRITE OUTPUT BED FILE ---
    list.beds[[feat]] <- file.path(dir.data, paste("kwtest_", feat, ".bed", sep=""))
    write.table(bed, list.beds[[feat]], sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
}




#############################################################################################################################
### GET FILES ---
files.bigwig <- sapply(subtypes, function(x) { paste(dir.bigwig, "/", sprintf("%s", x), ".bw", sep="") })
files.bigwig <- paste(files.bigwig, collapse=" ")

files.bed <- paste(unlist(list.beds), collapse=" ")


### COMPUTE MATRIX ---
file.matrix <- file.path(dir.data, "kwtest.tsv.gz")
cmd1 <- getComuteMatrix(files.bigwig, files.bed, file.matrix, sampleids=paste(subtypes, collapse=" "), n_cores=50)
system(cmd1)


### PLOT HEATMAP ---
file.plot <- file.path(dir.reproduce_fig, "figure_3_c.pdf")
cmd2 <- getPlotHeatmap(file.matrix=file.matrix, file.plot=file.plot)
system(cmd2)

