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

### DEFINE FILES ----
file.rds_atacseq_masterdata <- file.path(dir.reproduce_data, "wcdt_atacseq_masterdata.rds")
file.utility_functions <- file.path(dir.reproduce_scripts, "00_prepare_data/00_01_utility_functions.R")
file.annotation_gtf <- file.path(dir.wrk, "reference/annotation_data/gene_annotation/gencode_v28/gencode.v28.annotation.gtf.gz")

### LOAD FUNCTIONS ---
source(file.utility_functions)

### CHROMOSOMES ---
chromosomes <- paste("chr", c(1:22,"X","Y"), sep="")


###########################################################################################################################################
### FUNCTION: getFeatureRanges() ---
getFeatureRanges <- function(feature_ids){
    list.fids <- stringr::str_split(feature_ids, "_")
    chr <- unlist( lapply(list.fids, function(x) x[1]) )
    start <- as.numeric( unlist( lapply(list.fids, function(x) x[2]) ) )
    end <- as.numeric( unlist( lapply(list.fids, function(x) x[3]) ) )

    # CONVERT TO GENOMIC RANGES ---
    gr_peaks <- GenomicRanges::makeGRangesFromDataFrame(data.frame(chr=chr, start=start, end=end, feature_ids=feature_ids), keep.extra.columns=TRUE)

    return(gr_peaks)
}


### FUNCTION: annotPeaksbyRange() ---
annotPeaksbyRange <- function(gr_peaks, genecode.txdb){
    # ANNOTATE PEAKS ---
    annot_peaks <- ChIPseeker::annotatePeak(peak=gr_peaks, TxDb=genecode.txdb, verbose=TRUE)
    annot_peaks <- as.data.frame(annot_peaks)

    # ADD FEATURES ---
    annot_peaks$Feature <- ""
    annot_peaks$Feature[which(annot_peaks$annotation == "Promoter (<=1kb)")] <- "Promoter"
    annot_peaks$Feature[which(annot_peaks$annotation == "Promoter (1-2kb)")] <- "Promoter"
    annot_peaks$Feature[which(annot_peaks$annotation == "Promoter (2-3kb)")] <- "Promoter"
    annot_peaks$Feature[which(annot_peaks$annotation == "5' UTR")] <- "5' UTR"
    annot_peaks$Feature[which(annot_peaks$annotation == "3' UTR")] <- "3' UTR"

    exon <- annot_peaks$annotation[which(stringr::str_detect(annot_peaks$annotation, "Exon") == TRUE)]
    annot_peaks$Feature[which(annot_peaks$annotation %in% exon[which(str_detect(exon, "exon 1 of") == TRUE)])] <- "Exon"
    annot_peaks$Feature[which(annot_peaks$annotation %in% exon[which(str_detect(exon, "exon 1 of") == FALSE)])] <- "Exon"

    intn <- annot_peaks$annotation[which(stringr::str_detect(annot_peaks$annotation, "Intron") == TRUE)]
    annot_peaks$Feature[which(annot_peaks$annotation %in% intn[which(str_detect(intn, "intron 1 of") == TRUE)])] <- "Intron"
    annot_peaks$Feature[which(annot_peaks$annotation %in% intn[which(str_detect(intn, "intron 1 of") == FALSE)])] <- "Intron"

    annot_peaks$Feature[which(stringr::str_detect(annot_peaks$annotation, "Downstream") == TRUE)] <- "Downstream (<=300)"
    annot_peaks$Feature[which(annot_peaks$annotation == "Distal Intergenic")] <- "Distal Intergenic"

    return(annot_peaks)
}


### FUNCTION: addGenes() ---
addGenes <- function(annot, dat){
    pb <- txtProgressBar(min = 0, max = nrow(dat), style = 3)
    dat$Gene <- ""
    dat$GeneType <- ""
    n <- nrow(dat)

    # ADD GENE SYBBOLS ---
    for(i in 1:n){
        index <- which(annot$EnsemblTranscriptID == dat$transcriptId[i])
        if(length(index) != 0){
            dat$Gene[i] <- annot$Gene[index]
            dat$GeneType[i] <- annot$GeneType[index]

            # UPDATE PROGRESS BAR ---
            setTxtProgressBar(pb, i)
        }
    }

    return(dat)
}


###########################################################################################################################################
### TILE GENOME ---
bins <- binGenome(genome=BSgenome.Hsapiens.UCSC.hg38, chromosomes=chromosomes, tile.size=100)

### BINS TO FEATURE_IDS ---
df.bins <- as.data.frame(bins, stringsAsFactors=FALSE)
dt.bins <- data.table::as.data.table(df.bins)

### GET FEATURE IDS ---
mod_cols <- c("seqnames","start","end") 
myfun <- function(y) paste(y$seqnames, y$start, y$end, sep = "_")
feature_ids <- dt.bins[, myfun(.SD), .SDcols = mod_cols]

### GET RANGES ---
gr_bins <- getFeatureRanges(feature_ids)



###########################################################################################################################################
### LOAD GENE ANNOTATION DATABASE FROM GFF/GTF ---
genecode.txdb <- GenomicFeatures::makeTxDbFromGFF(file=file.annotation_gtf, format="gtf", organism="Homo sapiens", chrominfo=seqinfo(TxDb.Hsapiens.UCSC.hg38.knownGene))


### ANNOTATE FEATUREIDS ---
annot_bins <- annotPeaksbyRange(gr_peaks=gr_bins, genecode.txdb)
annot_bins <- subset(annot_bins, select=c("feature_ids","Feature","distanceToTSS","geneId","transcriptId"))
colnames(annot_bins) <- c("FeatureID","Feature","distanceToTSS","geneId","transcriptId")



###########################################################################################################################################
### GET ANNOTATION DATA ---
annot <- parseGTF(file.gtf=file.annotation_gtf, feature.type="transcript")

### ADD GENES ---
df <- addGenes(annot=annot, dat=annot_bins)




###########################################################################################################################################
### SAVE OBJECT TO RDATA FILE ---
file.rds_annot_bins <- file.path(dir.reproduce_data, "genomic_bins_feature_annotation.rds")
saveRDS(object=df, file=file.rds_annot_bins)


