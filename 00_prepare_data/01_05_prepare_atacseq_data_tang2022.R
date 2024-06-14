###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### DEFINE LIBRARIES ---
library("stringr")
library("data.table")
library("ChIPQC")
library("GenomicRanges")
library("DESeq2")
library("sva")

### DEFINE PATH ---
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/analysis/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

dir.process_tang <- file.path(dir.wrk, "external_data/tang_atacseq_2022/processed")
dir.process_wcdt <- file.path("/data1/datasets_1/human_prostate_WCDT/atac_solid/processed/2020_09_22/processed")

### DEFINE FILES ----
file.rds_atacseq_masterdata <- file.path(dir.reproduce_data, "wcdt_atacseq_masterdata.rds")
file.utility_functions <- file.path(dir.reproduce_scripts, "00_prepare_data/00_01_utility_functions.R")

### LOAD FUNCTIONS ---
source(file.utility_functions)


#####################################################################################################################

### LOAD ATAC-SEQ MASTER DATA ---
list.rds_atacseq_masterdata <- readRDS(file=file.rds_atacseq_masterdata)

### GET METADATA ---
metadata_tang2022_wcdt <- list.rds_atacseq_masterdata$metadata_tang2022_wcdt
sampleids_tang <- metadata_tang2022_wcdt$Sample_ID[which(metadata_tang2022_wcdt$Dataset == "TANG2022")]
sampleids_wcdt <- metadata_tang2022_wcdt$Sample_ID[which(metadata_tang2022_wcdt$Dataset == "WCDT")]

### LOAD ATAC-SEQ PEAKS ---
list.peaks_wcdt <- list.rds_atacseq_masterdata$atacseq_peaks
list.peaks_tang <- list.rds_atacseq_masterdata$atacseq_peaks_tang2022

### GET PEAK GENOMIC RANGES ---
list.peaks_wcdt_gr <- lapply(list.peaks_wcdt, ChIPQC:::GetGRanges, simple = TRUE)
list.peaks_tang_gr <- lapply(list.peaks_tang, ChIPQC:::GetGRanges, simple = TRUE)

### COMBINE PEAKS TO GET UNIQUE PEAKSET (COMBINED PEAKSET) ---
list.peaks_gr <- c(list.peaks_wcdt_gr, list.peaks_tang_gr)


#####################################################################################################################

### FUNCTION: getNonOverlapPeaks() ---
getNonOverlapPeaks_modified <- function(list.peaks_gr){
    list.peaks_gr <- GenomicRanges::GRangesList(list.peaks_gr)   
    reduced <- GenomicRanges::reduce(unlist(list.peaks_gr))
    consensusIDs <- paste0("consensus_", seq(1, length(reduced)))
    mcols(reduced) <- do.call(cbind, lapply(list.peaks_gr, function(x) (reduced %over% x) + 0))

    reducedConsensus <- reduced
    mcols(reducedConsensus) <- cbind(as.data.frame(mcols(reducedConsensus)), consensusIDs)
    consensusIDs <- paste0("consensus_", seq(1, length(reducedConsensus)))

    return(reducedConsensus)
}

### GET NON-OVERLAPPING PEAKS ---
consensusToCount <- getNonOverlapPeaks_modified(list.peaks_gr)


### 4. RETAIN COMBINED PEAKSET WITH OCCURANCE IN AT LEAST 2 SAMPLES ---
occurrences <- GenomicRanges::elementMetadata(consensusToCount) %>% 
                                    as.data.frame %>% 
                                    dplyr::select(-consensusIDs) %>% 
                                    rowSums

consensusToCount <- consensusToCount[occurrences >= 2,]


#####################################################################################################################


### FEATURE COUNTS: TANG2022 ----
dat_count_tang <- getATACseqFeatureCount(dir.process=dir.process_tang, consensusToCount=consensusToCount, sampleids=sampleids_tang, nthreads=40)

### FEATURE COUNTS: WCDT ----
dat_count_wcdt <- getATACseqFeatureCount(dir.process=dir.process_wcdt, consensusToCount=consensusToCount, sampleids=sampleids_wcdt, nthreads=40)


### MERGE DATA ---
dat_count <- merge(dat_count_wcdt, dat_count_tang, by="feature_id")
rownames(dat_count) <- dat_count$feature_id
dat_count$feature_id <- NULL


#####################################################################################################################


### PERFORM DESEQ2 (LIBRARY-SIZE) NORMALIZATION ---
dds <- DESeq2::DESeqDataSetFromMatrix(countData=dat_count, colData=subset(metadata, select=c("Sample_ID","Biopsy_Site")), design= ~ Biopsy_Site)
dds <- DESeq2::DESeq(dds)

### GET NORMALIZED FEATURE COUNTS ---
fcount_norm <- get.NormFeatureExpression(dds)

### QUANTILE NORMALIZATION ----
dQnorm <- get.normalizeQuantile(dat=fcount_norm)

#####################################################################################################################

##### COMBAT BATCH CORRECTION BETWEEN TANG2022 AND WCDT ATAC-SEQ DATASETS -------
### PREPARE MODEL MATRIX ---
des <- metadata_tang2022_wcdt
rownames(des) <- metadata_tang2022_wcdt$Sample_ID
des$Sample_ID <- NULL
des$Sample_Name <- NULL
model_matrix <- model.matrix(~as.factor(Biopsy_Site), data=des)

### RUN COMBAT ---
dat_combat_dqnorm <- sva::ComBat(dat=dQnorm, batch=des$Dataset, mod=model_matrix, par.prior=TRUE)

### RENAME COLNAMES ---
colnames(dat_combat_dqnorm) <- metadata_tang2022_wcdt$Sample_Name

#####################################################################################################################


### SAVE OBJECT TO RDATA FILE ---
file.rds_atacseq_read_counts_tang2022_wcdt_combined <- file.path(dir.reproduce_data, "atacseq_read_counts_tang2022_wcdt_combined.rds")
saveRDS(object=dat_combat_dqnorm, file=file.rds_atacseq_read_counts_tang2022_wcdt_combined)
