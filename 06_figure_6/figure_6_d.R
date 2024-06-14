###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### DEFINE LIBRARIES ---
library("stringr")

### DEFINE PATH ---
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/analysis/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

dir.bed <- file.path(dir.reproduce_data, "tf_footprints/mcrpc_subtypes/bed")
dir.bigwig <- file.path(dir.reproduce_data, "bigwig")

dir.output <- tempdir() # REPLACE PATH

### SUBTYPES ---
subtype <- "ARnNEn"
tfs <- c("myc","znf263")


##########################################################################################################################
### FUNCTION: getComuteMatrix() ---
getComuteMatrix <- function(file_bigwig, file_bed, file_matrix, sampleids, n_cores){
    # 'computeMatrix' FUNCTION FROM 'deepTools'
    # https://deeptools.readthedocs.io/en/latest/index.html

    cmd <- paste("computeMatrix reference-point",
                "--scoreFileName", file_bigwig,
                "--regionsFileName", file_bed,
                "--outFileName", file_matrix,
                "--referencePoint", "center",
                "--beforeRegionStartLength", 3000,
                "--afterRegionStartLength", 3000,
                "--binSize", 5,
                "--samplesLabel", sampleids,
                "--missingDataAsZero",                
                "--verbose",
                "-p", n_cores,  
                sep=" ")
    
    return(cmd)
}

### FUNCTION: getPlotHeatmap() ---
getPlotHeatmap <- function(file_matrix, file_plot, label_regions, label_samples){
    # 'plotHeatmap' FUNCTION FROM 'deepTools'
    # https://deeptools.readthedocs.io/en/latest/index.html

    cmd <- paste("plotHeatmap",
                "--matrixFile", file_matrix,
                "--outFileName", file_plot,
                "--sortRegions", "descend",
                "--sortUsing", "mean",
                "--sortUsingSamples", "1 2",
                "--colorMap", "jet",
                "--interpolationMethod", "gaussian",
                "--whatToShow", "\"heatmap and colorbar\"",
                "--refPointLabel", "\"Center\"",               
                "--regionsLabel", label_regions,
                "--samplesLabel", label_samples, 
                "--legendLocation", "\"best\"",
                "--xAxisLabel", "\"Footprints\"",
                "--heatmapHeight 4",
                "--heatmapWidth 2",                
                "--dpi", 600,
                "--verbose",
                sep=" ")
    
    return(cmd)
}

##########################################################################################################################
### GET FILES ---
files.bw <- sapply(tfs, function(x) { paste(dir.bigwig, "/tobias_footprint_", sprintf("%s_", x), subtype, ".bw", sep="") })
files.bed <- sapply(tfs, function(x) { paste(dir.bed, "/footprint_extended_", sprintf("%s_", x), subtype, ".bed", sep="") })



### FOOTPRINTS: ZNF263 ----
file_bigwig <- paste(files.bw, collapse=" ")
file_bed_1 <- files.bed["znf263"]
file_matrix_1 <- tempfile(pattern="footprint_matrix_znf263_regions", tmpdir=tempdir(), fileext=".tsv.gz")
file_plot_1 <- file.path(dir.reproduce_fig, "figure_6_d_1.pdf")

# COMPUTE MATRIX ---
cmd11 <- getComuteMatrix(file_bigwig=file_bigwig, 
                        file_bed=file_bed_1, 
                        file_matrix=file_matrix_1, 
                        sampleids=paste(c("myc","znf263"), collapse=" "), 
                        n_cores=55)
system(cmd11)

# PLOT HEATMAP ---
cmd12 <- getPlotHeatmap(file_matrix=file_matrix_1, 
                        file_plot=file_plot_1, 
                        label_regions="ZNF263_footprints", 
                        label_samples=paste("MYC","ZNF263", collapse=" "))
system(cmd12)





### FOOTPRINTS: MYC ----
file_bigwig <- paste(files.bw, collapse=" ")
file_bed_2 <- files.bed["myc"]
file_matrix_2 <- tempfile(pattern="footprint_matrix_myc_regions", tmpdir=tempdir(), fileext=".tsv.gz")
file_plot_2 <- file.path(dir.reproduce_fig, "figure_6_d_2.pdf")

# COMPUTE MATRIX ---
cmd21 <- getComuteMatrix(file_bigwig=file_bigwig, 
                        file_bed=file_bed_2, 
                        file_matrix=file_matrix_2, 
                        sampleids=paste(c("myc","znf263"), collapse=" "), 
                        n_cores=55)
system(cmd21)

# PLOT HEATMAP ---
cmd22 <- getPlotHeatmap(file_matrix=file_matrix_2, 
                        file_plot=file_plot_2, 
                        label_regions="MYC_footprints", 
                        label_samples=paste("MYC","ZNF263", collapse=" "))
system(cmd22)






### MERGE PLOT PDFs ---------------------------------------------------
list.files_plot <- list(file_plot_1, file_plot_2)
file.merge_plot <- file.path(dir.reproduce_fig, "figure_6_d.pdf")
cmd_qpdf <- paste("qpdf --empty --pages",
                paste(unlist(list.files_plot), collapse=" "),
                "--", file.merge_plot,
                sep=" ")

system(cmd_qpdf)

cat("\n", "FILE GENERATED:", file.merge_plot, "\n", sep=" ")
