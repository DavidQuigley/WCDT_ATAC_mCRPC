###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### DEFINE LIBRARIES ---
library("stringr")
library("data.table")
library("ggplot2")
library("gridExtra")
library("rGREAT")

### DEFINE PATH ---
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/analysis/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

### DEFINE FILES ---
file.rds_mat <- file.path(dir.reproduce_data, "atacseq_peaks_overlap_count_by_phenotype.rds")
file.rds_fannot <- file.path(dir.reproduce_data, "genomic_bins_feature_annotation.rds")
file.rds_atacseq_chromatin_variants <- file.path(dir.reproduce_data, "atacseq_chromatin_variants_benign_localizedpca_mcrpc_adeno_nepc.rds")


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

### FUNCTION: processEnrichment() ---
processEnrichment <- function(dat, pval){
    dat <- subset(dat, select=c("id","fold_enrichment","p_value","p_adjust"))
    dat <- dat[which(dat$p_adjust <= pval),]
    dat <- dat[order(dat$p_value, decreasing=FALSE),]
    dat$nlogp <- -log10(dat$p_adjust)

    return(dat)
}

### FUNCTION: plotBar() ----
plotBar <- function(dat, analysis_id){
    # PREPARE DATA ---
    dat <- dat[order(dat$nlogp, decreasing=FALSE),]
    dat$id <- factor(dat$id, levels=dat$id)

    # PLOT ---
    p <- ggplot(dat, aes(x=nlogp, y=id)) +
            geom_bar(stat="identity", fill="#006094", width=0.8) +
            coord_cartesian(xlim=c(0,8)) +
            scale_x_continuous(breaks=seq(0,8,by=2)) +
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
                legend.text = element_text(size = 5, color="#000000"),
                legend.title = element_blank(),
                legend.key.size = unit(0.3, "cm"),
                legend.position = "none") +
            ylab("") +
            xlab("-log10(fdr)") + 
            ggtitle(analysis_id) 

    return(p)
}

##############################################################################################################
### LOAD ATAC-SEQ PEAK OVERLAP BY PHENOTYPE ---
fannot <- readRDS(file=file.rds_fannot)


### LOAD ATAC-SEQ PEAK OVERLAP BY PHENOTYPE ---
mat <- readRDS(file=file.rds_mat)
mat <- mat[which(rowSums(mat) != 0),]
mat <- as.data.frame(mat)
mat <- cbind(FeatureID=rownames(mat), mat)
rownames(mat) <- NULL

### ADD KEY ---
mat$key <- apply(mat, 1, function(x) paste(as.numeric(x[2]), as.numeric(x[3]), as.numeric(x[4]), as.numeric(x[5]), sep=""))

#> sort(table(mat$key))
#   1001    1101    1011    1100    1010    1000    0101    1110    0011    0001
#    553    1890    8440    9981   10863   13530   43323   56531  166359  205531
#   1111    0111    0110    0100    0010
# 255654  334349  501657  893011 1155998

##############################################################################################################

### EXTRACT BY COMMON FEATURE IDS AND ORDER DATA ---
featuresids_common <- intersect(mat$FeatureID, fannot$FeatureID)
fannot <- subset(fannot, fannot$FeatureID %in% featuresids_common)
mat <- subset(mat, mat$FeatureID %in% featuresids_common)

fannot <- fannot[match(featuresids_common, fannot$FeatureID),]
mat <- mat[match(featuresids_common, mat$FeatureID),]



### MERGE ANNOTATION AND MAT ---
dat <- merge(mat, fannot, by="FeatureID")
gr_dat <- getGRangesObj(df=dat)




##############################################################################################################
### LOAD CHROMATIN VARIANTS (DIFF. PEAKS BETWEEN PHENOTYPES) ---
list.diff_peaks <- readRDS(file=file.rds_atacseq_chromatin_variants)
dat.diff_peaks <- list.diff_peaks$atacseq_diffpeaks_up_open
gr_diffPeaks <- getGRangesObj(df=dat.diff_peaks)


### FILTER OUT ALL PEAKS THAT HAVE DIFFERENTIAL ACCESSIBILITY BETWEEN STAGES ---
gr_unique <- IRanges::subsetByOverlaps(x=gr_dat, ranges=gr_diffPeaks, type="any", invert=TRUE)
dat_unique <- as.data.frame(gr_unique)
dat_unique <- subset(dat_unique, select=colnames(dat))

#############################################################################################
### [PCA + mCRPCAdeno + NEPC] i.e. 0111 
dat_pca <- dat_unique[which(dat_unique$key == "0111"),]
dat_pca <- subset(dat_pca, dat_pca$GeneType == "protein_coding")
dat_pca_promoter <- subset(dat_pca, dat_pca$Feature == "Promoter")
dat_pca_promoter_dtss0 <- subset(dat_pca_promoter, dat_pca_promoter$distanceToTSS == 0)

### LOAD PEAK BED FILE ---
gr_malignant <- getGRangesObj(df=dat_pca_promoter_dtss0)






##############################################################################################################
### SET SEED ---
set.seed(12345)

### GREAT ENRICHMENT ---
obj_reactome <- rGREAT::great(gr=gr_malignant, gene_sets="msigdb:C2:CP:REACTOME", tss_source="txdb:hg38", 
                            biomart_dataset = NULL, min_gene_set_size = 5, 
                            mode = "basalPlusExt", basal_upstream = 5000, basal_downstream = 1000, extension = 1000000, 
                            extended_tss = NULL, background = NULL, exclude = "gap", cores = 50, verbose = TRUE)

obj_hallmark <- rGREAT::great(gr=gr_malignant, gene_sets="MSigDB:H", tss_source="txdb:hg38", 
                            biomart_dataset = NULL, min_gene_set_size = 5, 
                            mode = "basalPlusExt", basal_upstream = 5000, basal_downstream = 1000, extension = 1000000, 
                            extended_tss = NULL, background = NULL, exclude = "gap", cores = 50, verbose = TRUE)



### PARSE REACTOME ---
dat_reactome <- obj_reactome@table
df_reactome <- processEnrichment(dat=dat_reactome, pval=0.05)
df_reactome$nlogp[1] <- df_reactome$nlogp[3] + 0.2
df_reactome$nlogp[2] <- df_reactome$nlogp[3] + 0.1

df_reactome <- df_reactome[1:40,]
del.index <- c(2,3,8,9,10,15,17,18,19,20,25,26,29,30,31,33,34,35,36,38,39,40)
df_reactome <- df_reactome[-del.index,]


### PARSE HALLMARK ---
dat_hallmark <- obj_hallmark@table
df_hallmark <- processEnrichment(dat=dat_hallmark, pval=0.05)


### PLOT  ---
p_reactome <- plotBar(dat=df_reactome, analysis_id="Reactome")
p_hallmark <- plotBar(dat=df_hallmark, analysis_id="Hallmark")



################################################################################
### GENERATE PLOT ---
file.plot <- file.path(dir.reproduce_fig, "supplementary_figure_06_a.pdf")
pdf(file.plot, width=4, height=2)
    grid.arrange(p_reactome, ncol=1, nrow=1)
dev.off()

file.plot <- file.path(dir.reproduce_fig, "supplementary_figure_06_b.pdf")
pdf(file.plot, width=2.5, height=1.5)
    grid.arrange(p_hallmark, ncol=1, nrow=1)
dev.off()
