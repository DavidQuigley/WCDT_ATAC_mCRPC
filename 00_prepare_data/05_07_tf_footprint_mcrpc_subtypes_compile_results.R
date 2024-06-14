###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### LOAD LIBRARIES ---
library("stringr")
library("data.table")
library("dplyr")

### DEFINE PATH ---
### DEFINE PATH ---
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/analysis/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

#dir.footprints <- file.path(dir.wrk, "analysis/12_tf_footprint/subtypes_AR/tobias")
dir.footprints <- file.path(dir.reproduce_data, "tf_footprints/mcrpc_subtypes")
dir.bam <- file.path(dir.reproduce_data, "tf_footprints/mcrpc_subtypes/bam")
dir.bed <- file.path(dir.reproduce_data, "tf_footprints/mcrpc_subtypes/bed")

### DEFINE FILES ---
file.motifdb <- file.path(dir.reproduce, "reference/JASPAR2018_motif_gene_conversion.tsv")

### SUBTYPES ---
subtypes <- c("ARpNEn","ARpNEp","ARlNEn","ARnNEp","ARnNEn")

### GET PAIRWISE ALL GROUP COMBINATIONS ---
dm_grp <- do.call(rbind.data.frame, combn(x=subtypes, m=2, simplify=FALSE)) %>%
                dplyr::mutate_all(as.character)
colnames(dm_grp) <- c("Group1","Group2")
rownames(dm_grp) <- NULL
dm_grp$id <- apply(dm_grp, 1, function(x) paste(x[1], x[2], sep="_"))
dm_grp$f_name <- apply(dm_grp, 1, function(x) paste("tobias_diff", paste(x[1], x[2], sep="_"), sep="_"))

### FUNCTION: ParseMotifDB() ---
ParseMotifDB <- function(file.motifdb){
    # LOAD MOTIFDB DATA ---
    motifdb <- read.delim(file.motifdb, header=TRUE, stringsAsFactors=FALSE)
    motifdb <- subset(motifdb, motifdb$species == "Homo sapiens")
    
    motifdb$motif_prefix <- stringr::str_replace_all(motifdb$motif_prefix, "::", "")
    motifdb$motif_prefix <- stringr::str_replace_all(motifdb$motif_prefix, "[(]", "")
    motifdb$motif_prefix <- stringr::str_replace_all(motifdb$motif_prefix, "[)]", "")

    motif_prefix <- stringr::str_replace_all(motifdb$motif_name, "::", "_")
    motif_prefix <- stringr::str_replace_all(motif_prefix, "[(]", "_")
    motif_prefix <- stringr::str_replace_all(motif_prefix, "[)]", "")
    motif_prefix <- stringr::str_replace_all(motif_prefix, "[.]", "_")

    motifdb$motif_prefix_2 <- paste(motifdb$motif_id, motif_prefix, sep="_")
    motifdb <- motifdb[order(motifdb$motif_id, decreasing=FALSE),]

    rownames(motifdb) <- NULL

    motifdb <- subset(motifdb, select=c("motif_id","motif_name","motif_prefix_2"))
    colnames(motifdb) <- c("motif_id","motif_name","motif_prefix")
    motifdb <- motifdb[!duplicated(motifdb),]

    return(motifdb)
}

### GET MOTIFS ---
motif_db <- ParseMotifDB(file.motifdb)



####################################################################################################################################################
### GET BINDETECT RESULTS FILES ---
files.bindetect_results <- sapply(dm_grp$f_name, function(x) { paste(dir.footprints, "/", sprintf("%s/", x), "bindetect_results.txt", sep="") })
names(files.bindetect_results) <- dm_grp$id

### PARSE BINDETECT RESULTS DATA ---
list.bindetect_results <- list()
for(k in 1:length(files.bindetect_results)){
    id <- names(files.bindetect_results)[k]

    # LOAD DATA ---
    dat <- data.table::fread(files.bindetect_results[id], data.table=FALSE, header=TRUE, stringsAsFactors=FALSE, nThread=1)

    # EXTRACT RESULTS FOR HUMAN TFs ONLY ---
    dat <- subset(dat, dat$motif_id %in% motif_db$motif_id)
    rownames(dat) <- NULL
    list.bindetect_results[[id]] <- dat
    
    cat("PROCESSED:", id, "\n", sep="\t")
}


### SAVE OBJECT TO RDATA FILE ---
file.rds_bindetect_results <- file.path(dir.footprints, "data_footprints/tobias_mcrpc_subtypes_bindetect_results.rds")
saveRDS(object=list.bindetect_results, file=file.rds_bindetect_results)



####################################################################################################################################################
### LOAD BINDETECT RESULTS ---
file.rds_bindetect_results <- file.path(dir.footprints, "data_footprints/tobias_mcrpc_subtypes_bindetect_results.rds")
list.bindetect_results <- readRDS(file=file.rds_bindetect_results)

### LOOP FOR EACH PAIR-WISE COMPARISON ---
list.footprints <- list()
for(k in 1:nrow(dm_grp)){
    id <- dm_grp$id[k]
    f_name <- dm_grp$f_name[k]
    subtype_1 <- dm_grp$Group1[k]
    subtype_2 <- dm_grp$Group2[k]

    cat("START:", id, "\n", sep="\t")

    # GET FOOTPRINTS OVERVIEW FILENAMES ---
    files.footprints_overview <- sapply(list.bindetect_results[[id]]$output_prefix, function(x) { paste(dir.footprints, "/", f_name, "/", sprintf("%s/%s", x, x), "_overview.txt", sep="") })
    names(files.footprints_overview) <- unlist(lapply(stringr::str_split(names(files.footprints_overview), "_"), function(x) x[1]))

    # LOAD FOOTPRINT OVERVIEW FILE ---
    list.dat <- list()
    for(i in 1:length(files.footprints_overview)){
        motif_id <- names(files.footprints_overview)[i]

        # LOAD DATA ---
        list.dat[[motif_id]] <- data.table::fread(files.footprints_overview[motif_id], data.table=FALSE, header=TRUE, stringsAsFactors=FALSE, nThread=1)

        cat("PROCESSED:", id, motif_id, "\n", sep=" ")
    }

    list.footprints[[id]] <- list.dat
    cat("PROCESSED:", id, "\n", sep="\t")
}



### SAVE OBJECT TO RDATA FILE ---
file.rds_tf_footprints_overview <- file.path(dir.footprints, "data_footprints/tobias_mcrpc_subtypes_tf_footprints_overview.rds")
saveRDS(object=list.footprints, file=file.rds_tf_footprints_overview)


