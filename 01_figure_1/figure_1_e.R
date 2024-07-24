###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### CONFIGURANTION ---
options("scipen"=10)

### DEFINE LIBRARIES ---
library("stringr")
library("data.table")
library("GenomicRanges")
library("IRanges")
library("reshape2")
library("ggplot2")
library("gridExtra")

### DEFINE PATH ---
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

dir.output <- tempdir() # REPLACE PATH
dir.bigwig <- file.path(dir.wrk, "secondary_data/wcdt_atacseq/bigwig_merge")  # REPLACE PATH

### DEFINE FILES ----
file.rds_atacseq_chromatin_variants <- file.path(dir.reproduce_data, "atacseq_chromatin_variants_benign_localizedpca_mcrpc_adeno_nepc.rds")
file.rds_atacseq_read_counts_combined <- file.path(dir.reproduce_data, "atacseq_read_counts_combined.rds")
file.rds_chipseq_bed <- file.path(dir.reproduce_data, "chipseq_bed.rds")


### LOAD ATAC-SEQ READ COUNTS DATA ---
list.rds_atacseq_read_counts_combined <- readRDS(file=file.rds_atacseq_read_counts_combined)

### GET FEATURE ANNOTATION DATA ---
dat.annot <- list.rds_atacseq_read_counts_combined$feature_annotation
dat.annot <- subset(dat.annot, select=c("FeatureID","Feature","seqnames","start","end"))

### LOAD ATAC-SEQ CHROMATIN VARIANTS DATA ---
list.rds_atacseq_chromatin_variants <- readRDS(file=file.rds_atacseq_chromatin_variants)

### GET CLUSTERED FEATURES ---
dat_atacseq <- list.rds_atacseq_chromatin_variants$atacseq_diffpeaks_cluster_features
dat_atacseq <- subset(dat_atacseq, dat_atacseq$ClusterName == "mcrpc_adeno")



### MERGE DATA ---
dat_atacseq <- merge(dat_atacseq, dat.annot, by="FeatureID")
dat_atacseq <- subset(dat_atacseq, dat_atacseq$Feature %in% c("Promoter","Intron","Distal Intergenic"))
dat_atacseq$Feature[which(dat_atacseq$Feature == "Distal Intergenic")] <- "Intergenic"

### GET GRANGES OF ATAC-SEQ CHROMATIN VARIANTS IN mCRPC-Adeno ---
gr_atacseq <- GenomicRanges::makeGRangesFromDataFrame(dat_atacseq, keep.extra.columns=TRUE, ignore.strand=TRUE)



### LOAD CHIP-SEQ DATA ---
list.rds_chipseq_bed <- readRDS(file=file.rds_chipseq_bed)
gr_chipseq <- list.rds_chipseq_bed$pomerantz_mcrpc$H3K27ac


### GET OVERLAP ---
gr_overlap <- IRanges::subsetByOverlaps(x=gr_atacseq, ranges=gr_chipseq, type="any")
df_overlap <- as.data.frame(gr_overlap)


### SAGGREGATE BY FEATURES ---
df_overlap_promoter <- subset(df_overlap, df_overlap$Feature == "Promoter")
df_overlap_intron <- subset(df_overlap, df_overlap$Feature == "Intron")
df_overlap_intergenic <- subset(df_overlap, df_overlap$Feature == "Intergenic")

### GET BED ---
features <- c("promoter","intron","intergenic")
list.bed <- list()
list.bed[["promoter"]]  <- df_overlap_promoter[,1:3]
list.bed[["intron"]]  <- df_overlap_intron[,1:3]
list.bed[["intergenic"]]  <- df_overlap_intergenic[,1:3]


### OUTPUT BED FILES ---
list.files_bed <- list()
list.files_matrix <- list()
for(feature in features){
    list.files_matrix[[feature]] <- tempfile(pattern=paste("diffpeaks_h3k27ac", feature, sep="_"), tmpdir=tempdir(), fileext=".tsv.gz")
    list.files_bed[[feature]] <- tempfile(pattern=paste("diffpeaks_h3k27ac", feature, sep="_"), tmpdir=tempdir(), fileext=".bed")
    write.table(list.bed[[feature]], list.files_bed[[feature]], sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
}



### FUNCTION: getComuteMatrix() ---
getComuteMatrix <- function(files.bigwig, file.bed, file.matrix, sampleids, n_cores=10){
    # 'computeMatrix' FUNCTION FROM 'deepTools'
    # https://deeptools.readthedocs.io/en/latest/index.html

    cmd <- paste("computeMatrix reference-point",
                "--scoreFileName", files.bigwig,
                "--regionsFileName", file.bed,
                "--outFileName", file.matrix,
                "--referencePoint", "center",
                "--beforeRegionStartLength", 3000,
                "--afterRegionStartLength", 3000,
                "--binSize", 10,
                "--samplesLabel", sampleids,
                "--missingDataAsZero",                
                "--verbose",
                "-p", n_cores,  
                sep=" ")
    
    return(cmd)
}


### COMPUTE MATRIX: ATAC-SEQ ---
file.bigwig <- file.path(dir.bigwig, "mcrpc_adeno.bw")
for(feature in features){
    cmd <- getComuteMatrix(files.bigwig=file.bigwig, 
                            file.bed=list.files_bed[[feature]], 
                            file.matrix=list.files_matrix[[feature]], 
                            sampleids="mcrpc", 
                            n_cores=40)
    system(cmd)
}





##############################################################################################################################################

### FUNCTION: loadMatrix() ---
loadMatrix <- function(file.matrix){
    # LOAD MATRIX ---
    dat <- data.table::fread(file=file.matrix, sep="\t", header=FALSE, nThread=1, skip=1, data.table=FALSE, verbose=FALSE)
    dat$V4  <- stringr::str_replace_all(stringr::str_replace_all(dat$V4, ":", "_"), "-", "_")
    rownames(dat) <- dat$V4
    dat <- dat[,-c(1:6)]

    # ADD COLNAMES ---
    bins <- c(seq(-3000, -10, by=10), seq(10,3000, by=10))
    colnames(dat) <- bins

    return(dat)
}

### FUNCTION: getMeanProfile() ---
getMeanProfile <- function(dat){
    # RESHAPE DATA ---
    df <- reshape2::melt(as.matrix(dat))
    colnames(df) <- c("FeatureID","Position","Value")

    # GET MEAN PROFILE ---
    dm <- df %>% 
            dplyr::group_by(Position) %>%
            dplyr::summarize(Mean=mean(Value))
    
    return(dm)
}

### FUNCTION: getPlot() ---
getPlot <- function(dm, title.main){
    # PLOT ---
    p <- ggplot(data=dm, aes(x=Position, y=Mean)) +       
                geom_line(aes(color=Feature), size=0.5, alpha=0.8) +
                scale_x_continuous(breaks=c(-3000,-2000,-1000,0,1000,2000,3000), labels=c("-3","-2","-1","C","1","2","3")) +
                scale_color_manual(values=c("#6a3d9a","#b15928","#e31a1c")) +
                #facet_grid(Phenotype~Gene, scales="free_y", space="fixed", shrink=TRUE) +
                #scale_y_continuous(breaks=seq(0,1.25, by=0.25)) + #ATAC-SEQ
                #coord_cartesian(ylim=c(0,1.25)) + #ATAC-SEQ
                theme(
                    axis.text.x = element_text(size = 5, color="#000000"),
                    axis.text.y = element_text(size = 5, color="#000000"),
                    axis.title = element_text(size = 5, color="#000000"),
                    plot.title = element_text(size = 5, color="#000000", hjust=0.5),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    axis.ticks = element_line(size=0.2, color="#000000"),	
                    strip.text = element_text(size=5, color="#000000"),
                    strip.background = element_rect(fill="#FFFFFF", color="#FFFFFF"),
                    panel.background = element_rect(fill="#FFFFFF", color="#000000"),
                    legend.text = element_text(size = 6, color="#000000"),
                    legend.title = element_blank(),
                    legend.key.size = unit(0.3, "cm"),
                    legend.position = "none") +
                ylab("ATAC-seq Profile") +
                xlab("Distance From Peak Center") + 
                ggtitle(title.main)
    
    return(p)
}


### LOAD MATRIX ---
list.dat <- list()
for(feature in features){
    dat <- loadMatrix(file.matrix=list.files_matrix[[feature]])
    dm <- getMeanProfile(dat)
    dm$Feature <- feature
    list.dat[[feature]] <- dm
}


### AGGREGATE DATA ---
dm <- do.call(rbind.data.frame, list.dat)
rownames(dm) <- NULL

#### FACTORIZE ---
dm$Feature <- factor(dm$Feature, levels=c("promoter","intron","intergenic"))

### PLOT ---
p <- getPlot(dm=dm, title.main="h3k27ac overlap atacseq diffpeaks")

# WRITE PLOT ---
file.plot <- file.path(dir.reproduce_fig, "figure_1_e.pdf")
pdf(file.plot, height=2, width=2)
    grid.arrange(p, nrow=1, ncol=1) 
dev.off()


