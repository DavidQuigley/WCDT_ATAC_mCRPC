###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### DEFINE LIBRARIES ---
library("stringr")
library("data.table")
library("dplyr")
library("tibble")
library("ggalluvial")
library("ggplot2")
library("gridExtra")

### DEFINE PATH ---
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/analysis/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

### DEFINE FILES ----
file.rds_atacseq_masterdata <- file.path(dir.reproduce_data, "wcdt_atacseq_masterdata.rds")
file.utility_functions <- file.path(dir.reproduce_scripts, "00_prepare_data/00_01_utility_functions.R")

### LOAD FUNCTIONS ---
source(file.utility_functions)

### CHROMOSOMES ---
chromosomes <- paste("chr", c(1:22,"X","Y"), sep="")


### FUNCTION: getBinsByPeaks() ---
getBinsByPeaks <- function(list.peaks, gr_bins, sampleid){
    # SGET PEAKS BY SAMPLEID ---
    df <- list.peaks[[sampleid]]

    # CONVERT TO GENOMIC RANGES ---
    gr_peaks <- GenomicRanges::makeGRangesFromDataFrame(df, keep.extra.columns=FALSE)
    gr_peaks <- GenomicRanges::sort(gr_peaks)
    gr_peaks <- IRanges::unique(gr_peaks)

    # SUBSET OVERLAPPING BIN REGIONS ---
    gr_bins_peaks <- IRanges::subsetByOverlaps(gr_bins, gr_peaks, ignore.strand=TRUE)

    return(gr_bins_peaks)
}

### FUNCTION: assignPeaks2Bin() ---
assignPeaks2Bin <- function(bins, feature_ids, mod_cols, sampleids, list.peaks){
    # FEATURE PEAK COUNTS MATRIX ---
    mat <- matrix(0, nrow=length(feature_ids), ncol=length(sampleids), dimnames=list(feature_ids, sampleids))

    # ASSIGN PEAKS TO BINS  ---
    for(sampleid in sampleids){
        cat("START:", sampleid, "\n", sep=" ")

        # INTERSECT BINS AND PEAKS ---
        gr_bin_peaks <- getBinsByPeaks(list.peaks, gr_bins=bins, sampleid)

        # EXTRACT FEATURE IDS ---
        df_bin_peaks <- as.data.frame(gr_bin_peaks)
        dt_bin_peaks <- data.table::as.data.table(df_bin_peaks)
        f_ids <- dt_bin_peaks[, myfun(.SD), .SDcols = mod_cols]

        # FLAG FEATURE IDS TO MATRIX ---
        mat[which(feature_ids %in% f_ids), sampleid] <- 1

        cat("PROCESSED:", sampleid, "\n", sep=" ")
    }

    return(mat)
}

#### getGroupFeatures() ---
getGroupFeatures <- function(mat){
    index <- which(rowSums(mat) != 0)
    f_ids <- rownames(mat)[index]
    return(f_ids)
}




### LOAD ATAC-SEQ MASTER DATA ---
list.rds_atacseq_masterdata <- readRDS(file=file.rds_atacseq_masterdata)

### GET SAMPLEIDS ---
sampleids_pome <- list.rds_atacseq_masterdata$metadata_combined$SAMPLE_ID[which((list.rds_atacseq_masterdata$metadata_combined$Dataset == "Pomerantz") & (list.rds_atacseq_masterdata$metadata_combined$Phenotype == "Benign") )]
sampleids_tcga <- list.rds_atacseq_masterdata$metadata_combined$SAMPLE_ID[which((list.rds_atacseq_masterdata$metadata_combined$Dataset == "TCGA") & (list.rds_atacseq_masterdata$metadata_combined$Phenotype == "PCa") )]
sampleids_mcrpc_adeno <- list.rds_atacseq_masterdata$metadata_mcrpc$Sample_Name[which(list.rds_atacseq_masterdata$metadata_mcrpc$mCRPC_Subtype != "AR-_NE+")]
sampleids_mcrpc_nepc <- list.rds_atacseq_masterdata$metadata_mcrpc$Sample_Name[which(list.rds_atacseq_masterdata$metadata_mcrpc$mCRPC_Subtype == "AR-_NE+")]


### GET PEAKS BY PHENOTYPES ---
list.atacseq_benign <- list.rds_atacseq_masterdata$atacseq_peaks_benign_pomerantz2020[sampleids_pome]
list.atacseq_pca <- list.rds_atacseq_masterdata$atacseq_peaks_localizedpca_tcgaprad[sampleids_tcga]
list.atacseq_mcrpc_adeno <- list.rds_atacseq_masterdata$atacseq_peaks[sampleids_mcrpc_adeno]
list.atacseq_mcrpc_nepc <- list.rds_atacseq_masterdata$atacseq_peaks[sampleids_mcrpc_nepc]




### TILE GENOME ---
bins <- binGenome(genome=BSgenome.Hsapiens.UCSC.hg38, chromosomes=chromosomes, tile.size=100)

### BINS TO FEATURE_IDS ---
df.bins <- as.data.frame(bins, stringsAsFactors=FALSE)
dt.bins <- data.table::as.data.table(df.bins)

### GET FEATURE IDS ---
mod_cols <- c("seqnames","start","end") 
myfun <- function(y) paste(y$seqnames, y$start, y$end, sep = "_")
feature_ids <- dt.bins[, myfun(.SD), .SDcols = mod_cols]




### GET PEAKS OVERLAP WITH BIN FEATURES ---
mat_benign <- assignPeaks2Bin(bins, feature_ids, mod_cols, sampleids=sampleids_pome, list.peaks=list.atacseq_benign)
mat_pca <- assignPeaks2Bin(bins, feature_ids, mod_cols, sampleids=sampleids_tcga, list.peaks=list.atacseq_pca)
mat_mcrpc_adeno <- assignPeaks2Bin(bins, feature_ids, mod_cols, sampleids=sampleids_mcrpc_adeno, list.peaks=list.atacseq_mcrpc_adeno)
mat_mcrpc_nepc <- assignPeaks2Bin(bins, feature_ids, mod_cols, sampleids=sampleids_mcrpc_nepc, list.peaks=list.atacseq_mcrpc_nepc)

### GET OVERLAP FEATURES ---
features_benign <- getGroupFeatures(mat=mat_benign)
features_pca <- getGroupFeatures(mat=mat_pca)
features_mcrpc_adeno <- getGroupFeatures(mat=mat_mcrpc_adeno)
features_mcrpc_nepc <- getGroupFeatures(mat=mat_mcrpc_nepc)

### CREATE FEATURE MATRIX BY PHENOTYPE---
phenotypes <- c("benign","localized_pca","mcrpc_adeno","mcrpc_nepc")
mat <- matrix(0, nrow=length(feature_ids), ncol=length(phenotypes), dimnames=list(feature_ids, phenotypes))

### FILL MATRIX ---
mat[which(rownames(mat) %in% features_benign), "benign"] <- 1
mat[which(rownames(mat) %in% features_pca), "localized_pca"] <- 1
mat[which(rownames(mat) %in% features_mcrpc_adeno), "mcrpc_adeno"] <- 1
mat[which(rownames(mat) %in% features_mcrpc_nepc), "mcrpc_nepc"] <- 1


### SAVE OBJECT TO RDATA FILE ---
file.rds_peak_overlap <- file.path(dir.reproduce_data, "atacseq_peaks_overlap_count_by_phenotype.rds")
saveRDS(object=mat, file=file.rds_peak_overlap)



### REMOVE EMPTY BINS ---
mat <- mat[which(rowSums(mat) != 0),]

### PREPARE DATA ---
mat <- as.data.frame(mat)
dat <- cbind(FeatureID=rownames(mat), mat)
rownames(dat) <- NULL


### PREPARE COUNT DATA ---
df <- dat %>% 
        dplyr::mutate_all(as.character) %>%
        #tibble::as_tibble(dat) %>%
        dplyr::count(benign,localized_pca,mcrpc_adeno,mcrpc_nepc, sort = FALSE, name = "Freq") %>%
        dplyr::mutate(benign=factor(benign, levels=c("0","1")), 
                    localized_pca=factor(localized_pca, levels=c("0","1")),
                    mcrpc_adeno=factor(mcrpc_adeno, levels=c("0","1")),
                    mcrpc_nepc=factor(mcrpc_nepc, levels=c("0","1")))


### FUNCTION: getAlluvialPlot() ---
getAlluvialPlot <- function(df){
    p <- ggplot(df, aes(y=Freq, axis1=benign, axis2=localized_pca, axis3=mcrpc_adeno, axis4=mcrpc_nepc)) +
            ggalluvial::geom_alluvium(aes(fill = mcrpc_adeno), width = 1/12) +
            ggalluvial::geom_stratum(width = 1/4, fill = "#FFFFFF", color = "#000000") +
            #geom_text(stat = "stratum", discern = TRUE, aes(label = after_stat(stratum))) +
            scale_x_discrete(limits = c("benign", "localized_pca", "mcrpc_adeno", "mcrpc_nepc"), expand = c(0.05, 0.05)) +
            scale_fill_manual(values = c("#377eb8","#e41a1c")) +
            theme(
			    axis.text.x = element_text(size = 5, color="#000000"),
			    axis.text.y = element_text(size = 5, color="#000000"),
			    axis.title = element_text(size = 5, color="#000000"),
			    plot.title = element_text(size = 5, color="#000000", hjust=0),
			    panel.grid.major = element_blank(),
			    panel.grid.minor = element_blank(),
			    axis.ticks = element_line(size=0.2, color="#000000"),
			    panel.background = element_rect(fill="#FFFFFF", color="#000000"),
			    legend.text = element_text(size = 5, color="#000000"),
			    legend.title = element_blank(),
			    legend.key.size = unit(0.2, "cm"),			
			    legend.position="none") + 
			ylab("") +
			xlab("") + 
            ggtitle("") 

    return(p)
}


### GENERATE ALLUVIAL PLOT ---
p <- getAlluvialPlot(df)

### WRITE PLOT ---
file.plot <- file.path(dir.reproduce_fig, "figure_1_b.pdf")
pdf(file.plot, width=2, height=2)  
    grid.arrange(p, nrow=1, ncol=1)
dev.off()
