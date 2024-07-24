###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### DEFINE LIBRARIES ---
library("stringr")
library("data.table")
library("IRanges")
library("GenomicRanges")
library("reshape2")
library("ggplot2")
library("gridExtra")

### DEFINE PATH ---
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

### DEFINE FILES ---
file.rds_mat <- file.path(dir.reproduce_data, "atacseq_peaks_overlap_count_by_phenotype.rds")
file.rds_atacseq_chromatin_variants <- file.path(dir.reproduce_data, "atacseq_chromatin_variants_benign_localizedpca_mcrpc_adeno_nepc.rds")
file.rds_great_enrichment_diffpeaks <- file.path(dir.reproduce_data, "great_enrichment_diffpeaks.rds") # OUTPUT RESULTS OF GREAT ENRICHMENT ANALYSIS

### DEFINE STAGES ---
keys <- c("1000","0100","0010","0001")
stages <- c("benign","pca","mcrpc","nepc")

##############################################################################################################
### FUNCTION: getGRangesObj() ---
getGRangesObj <- function(df){
    # GET GENOMIC COORDINATES ---
    df$chr <- unlist(lapply(stringr::str_split(df$FeatureID, "_"), function(x) x[1]))
    df$start <- as.numeric(unlist(lapply(stringr::str_split(df$FeatureID, "_"), function(x) x[2])))
    df$end <- as.numeric(unlist(lapply(stringr::str_split(df$FeatureID, "_"), function(x) x[3])))

    # GET GENOMIC RANGES OBJECT ---
    gr <- GenomicRanges::makeGRangesFromDataFrame(df, keep.extra.columns=TRUE)

    return(gr)
}

### FUNCTION: getBED() ---
getBED <- function(df){
    # GET GENOMIC COORDINATES ---
    df$chr <- unlist(lapply(stringr::str_split(df$FeatureID, "_"), function(x) x[1]))
    df$start <- as.numeric(unlist(lapply(stringr::str_split(df$FeatureID, "_"), function(x) x[2])))
    df$end <- as.numeric(unlist(lapply(stringr::str_split(df$FeatureID, "_"), function(x) x[3])))

    # GET GENOMIC RANGES OBJECT ---
    bed <- subset(df, select=c("chr","start","end"))
    rownames(bed) <- NULL

    return(bed)
}

### FUNCTION: getOverlappingPeaks() ---
getOverlappingPeaks <- function(diffPeaks, mat, key){
    # SUBSET DIFFPEAKS BY KEY ---
    diffPeaks_new <- diffPeaks[which(diffPeaks$key == key),]
    gr_diffPeaks <- getGRangesObj(df=diffPeaks_new)

    # SUBSET MAT BY KEY ---
    mat_new <- mat[which(mat$key == key),]
    gr_mat <- getGRangesObj(df=mat_new)

    # GET OVERLAPS ---
    gr_stage <- IRanges::subsetByOverlaps(x=gr_diffPeaks, ranges=gr_mat, type="any", invert=FALSE)
    dat_stage <- as.data.frame(gr_stage)
    dat_stage <- subset(dat_stage, select=colnames(diffPeaks))

    return(dat_stage)
}

### FUNCTION: parseEnrichment() ---
parseEnrichment <- function(df.enrich){
    # TRIM DATA --
    dat <- subset(df.enrich, select=c("ID","Desc","BinomRank","BinomP","BinomBonfP","BinomFdrQ","Genes"))

    # PARSE GO ID ---
    dat$ID <- stringr::str_replace_all(dat$ID, ":", "_")

    # SCORE ---
    dat$nlogBinomBonfP <- -log10(dat$BinomBonfP)
    dat$nlogBinomFdrQ <- -log10(dat$BinomFdrQ)

    # ORDER DATA ---
    dat <- dat[order(dat$nlogBinomFdrQ, decreasing=TRUE),]

    return(dat)
}

##############################################################################################################
### LOAD ATAC-SEQ PEAK OVERLAP BY PHENOTYPE ---
mat <- readRDS(file=file.rds_mat)
mat <- mat[which(rowSums(mat) != 0),]
mat <- as.data.frame(mat)
mat <- cbind(FeatureID=rownames(mat), mat)
rownames(mat) <- NULL
colnames(mat) <- c("FeatureID","benign","pca","mcrpc","nepc")

### ADD KEY ---
mat$key <- apply(mat, 1, function(x) paste(as.numeric(x[2]), as.numeric(x[3]), as.numeric(x[4]), as.numeric(x[5]), sep=""))


##############################################################################################################
### LOAD CHROMATIN VARIANTS (DIFF. PEAKS BETWEEN PHENOTYPES) ---
list.diff_peaks <- readRDS(file=file.rds_atacseq_chromatin_variants)
dat.diff_peaks <- list.diff_peaks$atacseq_diffpeaks_up_open

### SUBSET BY PROTEIN CODING GENES ---
dat.diff_peaks <- subset(dat.diff_peaks, dat.diff_peaks$GeneType %in% "protein_coding")
dat.diff_peaks$key <- apply(dat.diff_peaks, 1, function(x) paste(as.numeric(x[6]), as.numeric(x[7]), as.numeric(x[8]), as.numeric(x[9]), sep=""))


##############################################################################################################
list.bed <- list()
### LOOP FOE EACH STAGE ---
for(i in 1:length(stages)){
    stage <- stages[i]
    key <- keys[i]

    cat("START:", stage, "\n", sep=" ")
    
    # GET OVERLAP OF DIFFPEAKS AND THAT ARE EXCLUSIVELY PRESENT IN THE STAGE ---
    dat_stage <- getOverlappingPeaks(diffPeaks=dat.diff_peaks, mat=mat, key)

    # GET BED FILES ---
    bed_stage <- getBED(df=dat_stage)

    # STORE DATA ---
    list.bed[[stage]] <- bed_stage
    
    cat("PROCESSED:", stage, "\n", sep=" ")
}


### CONVERT TO GENOMIC RANGES ---
#list.gr_bed <- lapply(list.bed, function(x) GenomicRanges::makeGRangesFromDataFrame(x, keep.extra.columns=TRUE))






##############################################################################################################
# Note: Using the BED files generated above, 
# GREAT enrichment analysis was performed using the web-tool: http://great.stanford.edu/public/html/
##############################################################################################################

### LOAD GREAT ENRICHMENT ---
list.enrich <- readRDS(file=file.rds_great_enrichment_diffpeaks)

### PARSE ENRICHMENT RESULTS ---
list.enrich_parsed <- lapply(list.enrich, function(x) parseEnrichment(x))

### SELECT SIGNIFICANT TOP-25 PATHWAYS ---
list.enrich_pass <- lapply(list.enrich_parsed, function(x) x[which(x$BinomFdrQ <= 1e-10),])

### GET PATHWAY:GO-ID LOOPUP TABLE ---
pathway_annot <- do.call(rbind.data.frame, lapply(list.enrich_pass, function(x) x[1:25,1:2]))
pathway_annot <- na.omit(pathway_annot)
pathway_annot <- pathway_annot[!duplicated(pathway_annot),]
rownames(pathway_annot) <- NULL




### CREATE MATRIX ---
mat_enrich <- matrix(0, nrow=length(pathway_annot$ID), ncol=length(stages), dimnames=list(pathway_annot$ID, stages) )


### FILL MATRIX ---
for(j in 1:length(stages)){
    stage <- stages[j]

    # GET ENRICH DATA ---
    df <- list.enrich_pass[[stage]]

    # LOOP BY IDS ---
    for(i in 1:nrow(mat_enrich)){
        id <- rownames(mat_enrich)[i]
        index <- which(df$ID == id)

        if(length(index) != 0){
            mat_enrich[i,j] <- as.numeric( df$nlogBinomFdrQ[index] )
        }
    }

    cat("PROCESSED:", stage, "\n", sep=" ")  
}


### MERGE DATA ---
dat <- cbind(pathway_annot, mat_enrich)
rownames(dat) <- NULL
dat <- dat[order(dat$mcrpc, dat$pca, dat$nepc, decreasing=TRUE),]


### SELECTED GO:IDS ---
goids <- c("GO_0060520","GO_0071394","GO_0045945","GO_0045198","GO_0035089",
			"GO_0033148","GO_0006930","GO_0045113","GO_0030859","GO_0060768",
			"GO_0007010","GO_0097485","GO_0008202","GO_0000904","GO_0048667",
			"GO_0050653","GO_0008610","GO_0046398","GO_0030182","GO_0048666",
			"GO_0000902","GO_0048812","GO_0007409")

### EXTRACT BY GO:IDS ---
dat <- subset(dat, dat$ID %in% goids)
dat <- dat[,-which(colnames(dat) == "benign")]


### RESHAPE DATA ---
dm <- reshape2::melt(dat, id.vars=c("ID","Desc"), variable.name="stage", value.name="nlogBinomFdrQ")
dm$nlogBinomFdrQ[which(dm$nlogBinomFdrQ == 0)] <- NA

### FACTORIZE ---
dm$ID <- factor(dm$ID, rev(goids))
dm$stage <- factor(dm$stage, c("pca","mcrpc","nepc"))


### PLOT ---
p <- ggplot(dm, aes(x=stage, y=ID, size=nlogBinomFdrQ)) +
            geom_point(color="#006094") +
            scale_x_discrete(position="top") +
            scale_size_continuous(range = c(0, 2)) + # Adjust as required.
            theme(
                axis.text.x = element_text(size = 5, color="#000000"),
                axis.text.y = element_text(size = 5, color="#000000"),
                axis.title = element_text(size = 5, color="#000000"),
                plot.title = element_text(size = 5, color="#000000", hjust=0.5),
                panel.grid.major.y = element_line(size=0.1, color="#BDBDBD"),
                panel.grid.major.x = element_line(size=0.1, color="#BDBDBD"),
                panel.grid.minor = element_blank(),
                axis.ticks = element_line(size=0.2, color="#000000"), 
                strip.text = element_text(size=5, color="#000000"),
                strip.background = element_rect(fill="#FFFFFF", color="#FFFFFF"),
                panel.background = element_rect(fill="#FFFFFF", color="#000000"),
                legend.text = element_text(size = 5, color="#000000"),
                legend.title = element_blank(),
                legend.key.size = unit(0.3, "cm"),
                legend.position = "none") +
            ylab("") +
            xlab("") + 
            ggtitle("") 

# WRITE PLOT ---
file.plot <- file.path(dir.reproduce_fig, "figure_1_f.pdf")
pdf(file.plot, height=2, width=2)
    grid.arrange(p, nrow=1, ncol=1) 
dev.off()


