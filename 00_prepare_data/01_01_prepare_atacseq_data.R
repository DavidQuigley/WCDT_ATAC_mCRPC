###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### DEFINE LIBRARIES ---
library("stringr")
library("data.table")
library("dplyr")
library("tibble")

### DEFINE PATH ---
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

dir.pome <- file.path(dir.wrk, "external_data/pomerantz_atacseq_2020/processed")
dir.tcga <- file.path(dir.wrk, "external_data/tcgaprad_atacseq/processed")
dir.tang <- file.path(dir.wrk, "external_data/tang_atacseq_2022/processed")
dir.poppy <- file.path("/data1/projects/WCDT_WGS_2021/poppy")

### DEFINE FILES ---
file.metadata_mcrpc <- file.path(dir.metadata, "metadata_wcdt_atacseq_reproduce.tsv")
file.metadata_combined <- file.path(dir.metadata, "metadata_atacseq_combined_reproduce.tsv")
file.metadata_tang <- file.path(dir.wrk, "external_data/tang_atacseq_2022/metadata/metadata_tang_atacseq_2022.tsv")
file.metadata_tang_wcdt <- file.path(dir.wrk, "external_data/tang_atacseq_2022/metadata/metadata_tang2022_wcdt_combined.tsv")
file.wgs_status <- file.path(dir.metadata, "metadata_wcdt_wgs_status.tsv")
file.utility_functions <- file.path(dir.reproduce_scripts, "00_prepare_data/00_01_utility_functions.R")

file.atacseq_stats_pome <- file.path(dir.wrk, "external_data/pomerantz_atacseq_2020/analysis/01_atac_qc/data/pomerantz_atacseq_reads_peaks_frip_stats.tsv")
file.atacseq_stats_tcga <- file.path(dir.wrk, "external_data/tcgaprad_atacseq/analysis/01_atac_qc/data/tcga_prad_atacseq_reads_peaks_frip_stats.tsv")

file.rds_atacseq_read_counts_mcrpc <- file.path(dir.reproduce_data, "atacseq_read_counts_mcrpc.rds")

### LOAD FUNCTIONS ---
source(file.utility_functions)

### CHROMOSOMES ---
chromosomes <- paste("chr", c(1:22,"X","Y"), sep="")




####### METADATA: -----------------------------------------------------------------------------------------------------------------
### LOAD METADATA TABLE: mCRPC ---
metadata_mcrpc <- data.table::fread(file=file.metadata_mcrpc, sep="\t", header=TRUE, nThread=1,  data.table=FALSE, verbose=FALSE)

### LOAD METADATA TABLE: COMBINED DATASET ---
metadata_combined <- data.table::fread(file=file.metadata_combined, sep="\t", header=TRUE, nThread=1,  data.table=FALSE, verbose=FALSE)

### LOAD METADATA TABLE: TANG2022 AND WCDT COMBINED ---
metadata_tang <- data.table::fread(file=file.metadata_tang, sep="\t", header=TRUE, nThread=1,  data.table=FALSE, verbose=FALSE)
metadata_tang_wcdt <- data.table::fread(file=file.metadata_tang_wcdt, sep="\t", header=TRUE, nThread=1,  data.table=FALSE, verbose=FALSE)


####### COMPILE WCDT WGS DATA: -----------------------------------------------------------------------------------------------------------------
### GET POPPY FILES ---
sampleids_wgs <- as.character(na.omit(metadata_mcrpc$SampleID_WGS))
files.poppy_rds <- sapply(sampleids_wgs, function(x) { paste(dir.poppy, "/", sprintf("%s", x), "_poppy.Rdata", sep="") })

### LOAD POPPY RDATA ---
list.SO <- list()
for(sampleid in sampleids_wgs){
    if(file.exists(files.poppy_rds[sampleid])){
        # LOAD RDATA ---
        load(files.poppy_rds[sampleid])

        # STORE DATA ---
        list.SO[[sampleid]] <- somatic
        somatic <- NULL
    }

    cat("LOADED:", sampleid, "\n", sep=" ")
}

### COMPILE STATS ---
stats_wgs_mcrpc <- tibble::tibble(data.frame(SampleID=names(list.SO) )) %>%
                        dplyr::mutate_all(as.character) %>%
                        dplyr::mutate(Purity=as.numeric(unlist(lapply(list.SO, function(x) x$purity))),
                                DiploidProportion=as.numeric(unlist(lapply(list.SO, function(x) x$diploidProportion))),
                                Ploidy=as.numeric(unlist(lapply(list.SO, function(x) x$ploidy))),
                                WholeGenomeDuplication=toupper(as.character(unlist(lapply(list.SO, function(x) x$wholeGenomeDuplication)))),
                                MSstatus=as.character(unlist(lapply(list.SO, function(x) x$msStatus))),
                                TMBperMB=as.numeric(unlist(lapply(list.SO, function(x) x$tmbPerMb))),
                                TMBsv=as.numeric(unlist(lapply(list.SO, function(x) x$svTumorMutationalBurden))),
                                n_bnd=as.numeric(unlist(lapply(list.SO, function(x) x$n_bnd))),
                                n_del=as.numeric(unlist(lapply(list.SO, function(x) x$n_del))),
                                n_dup=as.numeric(unlist(lapply(list.SO, function(x) x$n_dup))),
                                n_ins=as.numeric(unlist(lapply(list.SO, function(x) x$n_ins))),
                                n_inv=as.numeric(unlist(lapply(list.SO, function(x) x$n_inv))) )



### GET mCRPC WGS STATUS ---
status_wgs_mcrpc <- data.table::fread(file=file.wgs_status, sep="\t", header=TRUE, nThread=1,  data.table=FALSE, verbose=FALSE)



####### ATAC-SEQ PEAKS -------------------------------------------------------------------------------------------------------------------
### WCDT mCRPC ---
list.atacseq.peaks <- getPeakBEDList(sampleids=metadata_mcrpc$Sample_Name, dir.process=dir.process, peak_type="narrow", peaks_shift=TRUE)

### WCDT mCRPC: UNCORRECTED ---
list.atacseq.peaks_uncorrected <- getPeakBEDList(sampleids=metadata_mcrpc$Sample_Name, dir.process=dir.process, peak_type="narrow", peaks_shift=FALSE)

### POMERANTZ ET AL. NATURE GENETICS 2020. BENIGN AND PCa---
sampleids_pome_benign <- metadata_combined$SAMPLE_ID[which((metadata_combined$Dataset == "Pomerantz") & (metadata_combined$Phenotype == "Benign") )]
sampleids_pome_localizedpca <- metadata_combined$SAMPLE_ID[which((metadata_combined$Dataset == "Pomerantz") & (metadata_combined$Phenotype == "PCa") )]

list.atacseq.benign_pomerantz2020 <- getPeakBEDList(sampleids=sampleids_pome_benign, dir.process=dir.pome, peak_type="narrow", peaks_shift=TRUE)
list.atacseq.localizedpca_pomerantz2020 <- getPeakBEDList(sampleids=sampleids_pome_localizedpca, dir.process=dir.pome, peak_type="narrow", peaks_shift=TRUE)


### CORCES ET AL. SCIENCE 2018. TCGA-PRAD LOCALIZED PCa ---
sampleid_tcga <- metadata_combined$SAMPLE_ID[which((metadata_combined$Dataset == "TCGA") & (metadata_combined$Phenotype == "PCa") )]
sampleid_tcga <- setdiff(sampleid_tcga, "TCGA-YL-A9WK-01A")
list.atacseq.localizedpca_tcgaprad <- getPeakBEDList(sampleids=sampleid_tcga, dir.process=dir.tcga, peak_type="narrow", peaks_shift=TRUE)


### TANG ET AL. SCIENCE 2022: mCRPC cell lines/organoids/PDX ----
tang2022_gsm_ids <- c("GSM5823615","GSM5823616","GSM5823628","GSM5823629","GSM5823630","GSM5823631","GSM5823632","GSM5823633")
list.atacseq.tang2022_1 <- getPeakBEDList(sampleids=metadata_tang$SAMPLE_ID[which(!(metadata_tang$GEO_ID %in% tang2022_gsm_ids))], dir.process=dir.tang, peak_type="narrow", peaks_shift=TRUE)
list.atacseq.tang2022_2 <- getPeakBEDList(sampleids=tang2022_gsm_ids, dir.process=dir.tang, peak_type="narrow", peaks_shift=TRUE)
list.atacseq.tang2022 <- c(list.atacseq.tang2022_1, list.atacseq.tang2022_2)



### WCDT CONSENSUS PEAKS ------------------------------------------------------------------------------------------------------------------------
### LOAD ATAC-SEQ READ COUNTS DATA ---
list.rds_atacseq_read_counts_mcrpc <- readRDS(file=file.rds_atacseq_read_counts_mcrpc)
dat_consensus <- list.rds_atacseq_read_counts_mcrpc$feature_annotation



####### ATAC-SEQ QC STATS -------------------------------------------------------------------------------------------------------------------
### POMERANTZ ET AL. 2020 ---
atacseq_stats_pomerantz2020 <- data.table::fread(file=file.atacseq_stats_pome, sep="\t", header=TRUE, nThread=1,  data.table=FALSE, verbose=FALSE)

### CORCES ET AL. 2018 (TCGA-PRAD) ---
atacseq_stats_tcgaprad <- data.table::fread(file=file.atacseq_stats_tcga, sep="\t", header=TRUE, nThread=1,  data.table=FALSE, verbose=FALSE)
atacseq_stats_tcgaprad$Phenotype <- "PCa"




###################################################################################################################################################
### ADD TO LIST ---
list.output <- list(metadata_mcrpc=metadata_mcrpc,
                        metadata_combined=metadata_combined,
                        metadata_tang2022=metadata_tang,
                        metadata_tang2022_wcdt=metadata_tang_wcdt,
                        stats_wgs_mcrpc=stats_wgs_mcrpc,
                        status_wgs_mcrpc=status_wgs_mcrpc,
                        atacseq_peaks=list.atacseq.peaks,
                        atacseq_peaks_consensus=dat_consensus,
                        atacseq_peaks_uncorrected=list.atacseq.peaks_uncorrected,
                        atacseq_peaks_benign_pomerantz2020=list.atacseq.benign_pomerantz2020,
                        atacseq_peaks_localizedpca_pomerantz2020=list.atacseq.localizedpca_pomerantz2020,
                        atacseq_peaks_localizedpca_tcgaprad=list.atacseq.localizedpca_tcgaprad,
                        atacseq_peaks_tang2022=list.atacseq.tang2022,
                        atacseq_stats_pomerantz2020=atacseq_stats_pomerantz2020,
                        atacseq_stats_tcgaprad=atacseq_stats_tcgaprad)


### SAVE OBJECT TO RDATA FILE ---
file.rds_atacseq_masterdata <- file.path(dir.reproduce_data, "wcdt_atacseq_masterdata.rds")
saveRDS(object=list.output, file=file.rds_atacseq_masterdata)
