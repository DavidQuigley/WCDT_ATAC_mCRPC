###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### LOAD LIBRARIES ---
library("stringr")
library("data.table")
library("dplyr")
library("reshape2")

### DEFINE PATH ---
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/analysis/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

dir.footprints <- file.path(dir.reproduce_data, "tf_footprints/mcrpc_subtypes")
dir.bam <- file.path(dir.reproduce_data, "tf_footprints/mcrpc_subtypes/bam")
dir.bed <- file.path(dir.reproduce_data, "tf_footprints/mcrpc_subtypes/bed")

### DEFINE FILES ---
file.rds_footprints <- file.path(dir.footprints, "data_footprints/expt_bound_pvalue_footprint_summary.rds")
file.rds_chipseq <- file.path(dir.reproduce_data, "chipseq_bed.rds")
file.utility_functions <- file.path(dir.reproduce_scripts, "00_prepare_data/00_01_utility_functions.R")

### LOAD FUNCTIONS ---
source(file.utility_functions)

### SUBTYPES ---
subtypes <- c("ARpNEn","ARpNEp","ARlNEn","ARnNEp","ARnNEn")

### PVALUES ---
pvals <- c(0.5, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001)
names(pvals) <- c("p5","p1","p05","p01","p005","p001","p0005","p0001")

### MOTIFS ---
motifs <- c("ExtendedSite_AR","MA0148.3_FOXA1","MA0901.1_HOXB13")



### GET PAIRWISE ALL GROUP COMBINATIONS ---
dm_grp <- do.call(rbind.data.frame, combn(x=subtypes, m=2, simplify=FALSE)) %>%
                dplyr::mutate_all(as.character)
colnames(dm_grp) <- c("Group1","Group2")
rownames(dm_grp) <- NULL
dm_grp$Analysis <- paste(dm_grp$Group1, dm_grp$Group2, sep="_")



### FUNCTION: prepareData() ---
prepareData <- function(dm_grp, list.dat, subtypes, pval, motif){
    # EXTRACT DATA ---
    list.df <- list()
    for(k in 1:nrow(dm_grp)){
        analysis <- dm_grp$Analysis[k]

        # LOAD DATA ---
        dat <- list.dat[[pval]][[analysis]][[motif]]

        # GET FEATUREIDS ---
        dat$FeatureID <- apply(dat, 1, function(x) paste(x[1],as.numeric(x[2]),as.numeric(x[3]), sep="_"))

        # GET SUBTYPE GROUP ---
        val_grp1 <- stringr::str_split(colnames(dat)[11], "_")[[1]][1]
        val_grp2 <- stringr::str_split(colnames(dat)[12], "_")[[1]][1]

        # EXTRACT DATA BY SUBTYPE ---
        df1 <- subset(dat, select=c("FeatureID", "peak_id", colnames(dat)[13]))
        df2 <- subset(dat, select=c("FeatureID", "peak_id", colnames(dat)[14]))

        colnames(df1) <- c("FeatureID","peak_id","Status")
        colnames(df2) <- c("FeatureID","peak_id","Status")

        # REMOVE SITES WITHOUT ATAC-SEQ PEAKS ---
        df1 <- df1[stringr::str_detect(df1$peak_id, val_grp1),]
        df2 <- df2[stringr::str_detect(df2$peak_id, val_grp2),]

        # REMOVE PEAKID COLUMN ---
        df1$peak_id <- NULL
        df2$peak_id <- NULL

        # ADD GROUP ---
        df1$Subtype <- val_grp1
        df2$Subtype <- val_grp2

        # ROW BIND DATA ---
        list.df[[analysis]] <- rbind(df1, df2)
    }

    # AGGREGATE DATA ---
    df_combn <- do.call(rbind.data.frame, list.df)
    rownames(df_combn) <- NULL

    # REMOVE DUPLICATE ---
    df_combn <- df_combn[!duplicated(df_combn),]

    # GET UNIQUE SET OF FEATURES ---
    feature_ids <- unique(df_combn$FeatureID)

    # LOOP BY SUBTYPE ---
    list.dm <- list()
    for(k in 1:length(subtypes)){
        subtype <- subtypes[k]
        df3 <- subset(df_combn, df_combn$Subtype == subtype)
        f_ids <- df3$FeatureID[which(df3$Status == 1)]

        df4 <- data.frame(motif=motif, pval=pval, Subtype=subtype, FeatureID=feature_ids, Status=0)
        df4$Status[which(df4$FeatureID %in% f_ids)] <- 1

        list.dm[[subtype]] <- df4
    }

    # AGGREGATE DATA ---
    #dm <- do.call(rbind.data.frame, list.dm)
    #rownames(dm) <- NULL

    return(list.dm)
}

### FUNCTION: getGenomeRanges() ---
getGenomeRanges <- function(list.df, motif, pval, subtype){
    # GET DATA ---
    df <- list.df[[motif]][[pval]][[subtype]]

    # GET COORDINATES ---
    df$chr <- unlist(lapply(stringr::str_split(df$FeatureID, "_"), function(x) x[1]))
    df$start <- as.numeric(unlist(lapply(stringr::str_split(df$FeatureID, "_"), function(x) x[2])))
    df$end <- as.numeric(unlist(lapply(stringr::str_split(df$FeatureID, "_"), function(x) x[3])))

    # CONVERT TO GRANGES OBJECT ---
    gr <- GenomicRanges::makeGRangesFromDataFrame(df, keep.extra.columns=TRUE, ignore.strand=TRUE)

    return(gr)
}

### FUNCTION: get_overlap() ---
get_overlap <- function(gr.atacseq, gr.chipseq){
    # GET OVERLAP ---
    gr.overlap_atacseq <- IRanges::subsetByOverlaps(x=gr.atacseq, ranges=gr.chipseq, type="any")

    return(gr.overlap_atacseq)
}    

### FUNCTION: getConfusionMatrix() ---
getConfusionMatrix <- function(gr.atacseq, gr.overlap){
    # GET OVERLAPPED FEATUREIDS ---
    feature_ids <- as.character(gr.overlap$FeatureID)

    # CONVERT TO DATAFRAME ---
    df <- as.data.frame(gr.atacseq)
    df <- subset(df, select=c("motif","pval","Subtype","FeatureID","Status"))

    # TAG OVERLAPPED FEATUREIDS ---
    df$ChipseqHits <- 0
    df$ChipseqHits[which(df$FeatureID %in% feature_ids)] <- 1

    # ADD TAGS ---
    df$tag <- apply(df, 1, function(x) paste(x[5], x[6], sep="_"))

    # GET CONTINGENCY TABLE ----
    tbl <- as.matrix(table(df$tag))

    # GET CONFUSION MATRIX ---
    tbl1 <- c(tbl["1_1",1], tbl["1_0",1], tbl["0_1",1], tbl["0_0",1])
    mat <- matrix(tbl1, nrow=2, ncol=2, byrow=TRUE, dimnames=list(c("1","0"),c("1","0")) )

    return(mat)
}    


### LOAD ATAC-SEQ FOOTPRINTS DATA ---
list.rds_footprints <- readRDS(file=file.rds_footprints)

### LOAD CHIP-SEQ DATA ---
list.rds_chipseq <- readRDS(file=file.rds_chipseq)
gr_chipseq_AR <- list.rds_chipseq$pomerantz_mcrpc$AR
gr_chipseq_FOXA1 <- list.rds_chipseq$pomerantz_mcrpc$FOXA1
gr_chipseq_HOXB13 <- list.rds_chipseq$pomerantz_mcrpc$HOXB13


### GET BOUND MOTIFS BY SUBTYPE ---
list.df <- list()
for(i in 1:length(motifs)){
    motif <- motifs[i]

    # LOOP BY ANALYSIS ---
    list.df_pval <- list()
    for(j in 1:length(pvals)){
        pval <- names(pvals)[j]
        list.df_pval[[pval]] <- prepareData(dm_grp, list.dat=list.rds_footprints, subtypes, pval, motif)
        cat("PROCESSED:", motif, pval, "\n", sep=" ")
    } 
    list.df[[motif]] <- list.df_pval
}




### GET GENOMIC RANGES ---
list.gr_atacseq <- list()
for(i in 1:length(motifs)){
    motif <- motifs[i]

    # LOOP BY PVAL ---
    list.gr_pval <- list()
    for(j in 1:length(pvals)){
        pval <- names(pvals)[j]

        # LOOP BY SUBTYPE ---
        list.gr_pval_subtype <- list()
        for(k in 1:length(subtypes)){
            subtype <- subtypes[k]

            list.gr_pval_subtype[[subtype]] <- getGenomeRanges(list.df, motif, pval, subtype)

            cat("PROCESSED:", motif, pval, subtype, "\n", sep=" ")
        }
        list.gr_pval[[pval]] <- list.gr_pval_subtype
    }
    list.gr_atacseq[[motif]] <- list.gr_pval    
}


### GET OVERLAP OF TF FOOTPRINT WITH CHIPSEQ PEAKS ---
list.overlap <- list()
for(i in 1:length(motifs)){
    motif <- motifs[i]

    # LOOP BY PVAL ---
    list.overlap_pval <- list()
    for(j in 1:length(pvals)){
        pval <- names(pvals)[j]

        # LOOP BY ANALYSIS ---
        list.overlap_pval_subtype <- list()
        for(k in 1:length(subtypes)){
            subtype <- subtypes[k]

            gr.atacseq <- list.gr_atacseq[[motif]][[pval]][[subtype]]
        
            # GET ATAC-SEQ AND CHIP-SEQ OVERLAP ---
            if(motif == "ExtendedSite_AR"){
                gr.overlap <- get_overlap(gr.atacseq, gr.chipseq=gr_chipseq_AR)
            }else if(motif == "MA0148.3_FOXA1"){
                gr.overlap <- get_overlap(gr.atacseq, gr.chipseq=gr_chipseq_FOXA1)        
            }else if(motif == "MA0901.1_HOXB13"){
                gr.overlap <- get_overlap(gr.atacseq, gr.chipseq=gr_chipseq_HOXB13)        
            }

            # PARSE OVERLAP DATA TO GET CONFUSION MATRIX ---
            list.overlap_pval_subtype[[subtype]] <- getConfusionMatrix(gr.atacseq, gr.overlap)
            
            cat("PROCESSED:", motif, pval, subtype, "\n", sep=" ")
        }
        list.overlap_pval[[pval]] <- list.overlap_pval_subtype
    }
    list.overlap[[motif]] <- list.overlap_pval
}


### SAVE OBJECT TO RDATA FILE ---
file.rds_confusion_matrix <- file.path(dir.footprints, "data_footprints/expt_footprint_chipseq_confusion_matrix.rds")
saveRDS(object=list.overlap, file=file.rds_confusion_matrix)




### COMPUTE PERFORMANCE MATRIX ---
list.performance <- list()
for(i in 1:length(motifs)){
    motif <- motifs[i]

    # LOOP BY PVAL ---
    list.performance_pval <- list()
    for(j in 1:length(pvals)){
        pval <- names(pvals)[j]

        # LOOP BY ANALYSIS ---
        list.performance_pval_subtype <- list()
        for(k in 1:length(subtypes)){    
            subtype <- subtypes[k]
            c_mat <- list.overlap[[motif]][[pval]][[subtype]]

            # COMPUTE PERFORMANCE ---
            list.performance_pval_subtype[[subtype]] <- computePerformance(mat=c_mat)
            
            cat("PROCESSED:", motif, pval, subtype, "\n", sep=" ")
        }
        list.performance_pval[[pval]] <- list.performance_pval_subtype
    }
    list.performance[[motif]] <- list.performance_pval
}




### COMPILE PERFORMANCE DATA ---
list.dm <- list()
ctr <- 1
for(i in 1:length(motifs)){
    motif <- motifs[i]

    # LOOP BY PVAL ---
    for(j in 1:length(pvals)){
        pval <- names(pvals)[j]

        # LOOP BY SUBTYPE ---
        for(k in 1:length(subtypes)){    
            subtype <- subtypes[k]

            list.pmetric <- list.performance[[motif]][[pval]][[subtype]]
            list.pmetric$metric <- NULL

            # COMPILE DATA ---
            df <- cbind(data.frame(motif=motif, pval=pval, subtype=subtype), t(as.data.frame(unlist(list.pmetric))))

            list.dm[[ctr]] <- df
            ctr <- ctr + 1
        }
    }
}

### AGGREGATE DATA ---
dm <- do.call(rbind.data.frame, list.dm)
rownames(dm) <- NULL

### WRITE OUTPUT ---
file.dm <- file.path(dir.footprints, "data_footprints/expt_footprint_chipseq_prediction_performance_metrics.tsv")
write.table(dm, file.dm, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)


