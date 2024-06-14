###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### DEFINE LIBRARIES ---
library("stringr")
library("data.table")
library("Rsubread")
library("DESeq2")

### DEFINE PATH ---
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/analysis/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

dir.bam_wcdt_mcrpc <- file.path("/data1/datasets_1/human_prostate_WCDT/atac_solid/processed/2020_09_22/processed") #CHANGE TO AS REQUIRED

### DEFINE FILES ----
file.rds_atacseq_masterdata <- file.path(dir.reproduce_data, "wcdt_atacseq_masterdata.rds")
file.utility_functions <- file.path(dir.reproduce_scripts, "00_prepare_data/00_01_utility_functions.R")

### LOAD FUNCTIONS ---
source(file.utility_functions)




#####################################################################################################################
### LOAD ATAC-SEQ MASTER DATA ---
list.rds_atacseq_masterdata <- readRDS(file=file.rds_atacseq_masterdata)

### GET METADATA ---
metadata <- list.rds_atacseq_masterdata$metadata_mcrpc
sampleids <- metadata$Sample_ID


#####################################################################################################################
### GET ATAC-SEQ PEAKS: WCDT mCRPC ---
list.peaks <- list.rds_atacseq_masterdata$atacseq_peaks
names(list.peaks) <- metadata$Sample_ID


####################################################################################################
### GET NON-OVERLAPPING PEAKS: CONSENSUS PEAKS ---
# This function will provide a binary tag if the respective genomic region is covered in the sample or not
consensusToCount <- getNonOverlapPeaks(list.peaks, sampleids)


#> length(consensusToCount)
#[1] 682382

#> head(consensusToCount[,1:4])
#GRanges object with 6 ranges and 4 metadata columns:
#      seqnames      ranges strand | DTB.004.BL DTB.005.BL DTB.008.BL DTB.009.BL
#         <Rle>   <IRanges>  <Rle> |  <numeric>  <numeric>  <numeric>  <numeric>
#  [1]     chr1  9996-10678      * |          1          1          1          1
#  [2]     chr1 11139-11676      * |          1          1          1          1
#  [3]     chr1 12228-12382      * |          0          0          0          0
#  [4]     chr1 13408-13834      * |          0          0          0          0
#  [5]     chr1 13987-14240      * |          1          0          0          0
#  [6]     chr1 15501-15695      * |          0          0          0          0
#  -------
#  seqinfo: 24 sequences from an unspecified genome; no seqlengths


####################################################################################################
### RETAIN PEAKS THAT ARE DETEDCTED IN AT LEAST 2 SAMPLES ---
occurrences <- GenomicRanges::elementMetadata(consensusToCount) %>% 
                                    as.data.frame %>% 
                                    dplyr::select(-consensusIDs) %>% 
                                    rowSums

consensusToCount_atleast2samples <- consensusToCount[occurrences >= 2,]

#> length(consensusToCount_atleast2samples)
#[1] 355814




####################################################################################################
### BAM FILE PATHS: WCDT mCRPC ---
files.bam_wcdt_mcrpc <- sapply(metadata$Sample_Name, function(x) {paste(dir.bam_wcdt_mcrpc, "/", sprintf("%s/alignment/%s", x, x), ".bowtie.sorted.nodup.shifted.bam", sep="")}) #CHANGE TO AS REQUIRED
names(files.bam_wcdt_mcrpc) <- metadata$Sample_ID



### FUNCTION: getFeatureCount() ---
getFeatureCount <- function(files.bam, consensusToCount, sampleids, nthreads){
    # GET REGION TO COUNT ---
    regionsToCount <- data.frame(GeneID = paste(
                                            GenomicRanges::seqnames(consensusToCount), 
                                            GenomicRanges::start(consensusToCount), 
                                            GenomicRanges::end(consensusToCount), sep = "_"), 
                                Chr = GenomicRanges::seqnames(consensusToCount), 
                                Start = GenomicRanges::start(consensusToCount), 
                                End = GenomicRanges::end(consensusToCount), 
                                Strand = GenomicRanges::strand(consensusToCount))

    # COUNT READS
    fcResults <- Rsubread::featureCounts(
                                files=files.bam, 
                                annot.ext=regionsToCount, 
                                isPairedEnd=TRUE, 
                                countMultiMappingReads=FALSE, 
                                maxFragLength=100,
                                autosort=TRUE,
                                nthreads=nthreads,
                                verbose=FALSE)
  
    # COLLECT THE RAW COUNTS ---
    counts <- fcResults$counts
    colnames(counts) <- sampleids
    counts <- cbind.data.frame(feature_id=rownames(counts), counts)
    rownames(counts) <- NULL

    return(counts)
}

### RAW READ COUNTS ---
feature_counts_raw <- getFeatureCount(files.bam=files.bam_wcdt_mcrpc, consensusToCount=consensusToCount_atleast2samples, sampleids=metadata$Sample_ID, nthreads=50)

########################################################################################################
### NOTE: READ BELOW ---
# 'feature_counts_raw' is stored in the following RDS file:  file.path(dir.reproduce_data, "atacseq_read_counts_combined.rds")
########################################################################################################





########################################################################################################
##### DATA NORMALIZATION ###########
# 'feature_counts_raw': RAW READ COUNTS MATRIX
# 'colData': Data frame of sample_id and covariates (See DESeq2 manual: https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
# 'design': Formula or covariate to compare against (See DESeq2 manual: https://bioconductor.org/packages/release/bioc/html/DESeq2.html)

### PERFORM DESEQ2 (LIBRARY-SIZE) NORMALIZATION ---
dds <- DESeq2::DESeqDataSetFromMatrix(countData=feature_counts_raw, colData, design)
dds <- DESeq2::DESeq(dds)

### LOG TRANSFORM DESEQ dds OBJECT ---
rds <- DESeq2::vst(dds)
feature_counts_norm <- SummarizedExperiment::assay(rds)

########################################################################################################
### NOTE: READ BELOW ---
# 'feature_counts_norm' is stored in the following RDS file:  file.path(dir.reproduce_data, "atacseq_read_counts_combined.rds")
########################################################################################################

