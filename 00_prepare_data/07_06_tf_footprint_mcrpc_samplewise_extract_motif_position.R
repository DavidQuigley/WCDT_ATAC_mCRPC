###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### LOAD LIBRARIES ---
library("stringr")
library("data.table")

### DEFINE PATH ---
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

#dir.footprint_original <- file.path("/data1/projects/WCDT_atac_2020/analysis/12_tf_footprint/mcrpc_samplewise/tobias") # ORIGINAL DATA PRESENT HERE
dir.footprint <- file.path(dir.reproduce_data, "tf_footprints/mcrpc_samplewise")
dir.bed <- file.path(dir.footprint, "bed")

### DEFINE FILES ----
file.rds_atacseq_masterdata <- file.path(dir.reproduce_data, "wcdt_atacseq_masterdata.rds")
file.motifdb <- file.path(dir.reproduce, "reference/JASPAR2018_motif_gene_conversion.tsv")


### CHROMOSOMES ---
chromosome <- paste("chr", c(1:22,"X","Y"), sep="")



### FUNCTION: ParseMotifDB() ---
ParseMotifDB <- function(file.motifdb){
    # LOAD MOTIFDB DATA ---
    motifdb <- data.table::fread(file=file.motifdb, sep="\t", header=TRUE, nThread=1, data.table=FALSE, stringsAsFactors=FALSE, verbose=FALSE)
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

### FUNCTION: getFolderMotifdb() ---
getFolderMotifdb <- function(motif_db, dir.diff){
    # GET MOTIFS FOLDERS ---
    f_motifs <- list.dirs(dir.diff, recursive=FALSE, full.names=FALSE)
    f_ids <- unlist(lapply(stringr::str_split(f_motifs, "_"), function(x) x[1]))
    df_motifs <- data.frame(motif_folder=f_motifs, motif_id=f_ids) 
    df.motif_db <- merge(motif_db, df_motifs, by="motif_id")
    return(df.motif_db)
}

### FUNCTION: getMotifPosition() ---
getMotifPosition <- function(dir.footprint_input, sampleids, motif_prefix, motif_folder, bound_status){
    # GET FOOTPRINT FILES PATH ---
    files.footprint <- sapply(sampleids, function(x) { paste(dir.footprint_input, "/", sprintf("%s/", x), "bindetect/", motif_folder, "/", motif_folder, "_overview.txt", sep="") })

    # LOAD EACH FOOTPRINT FILE ---
    list.motif_sampleid <- list()
    for(k in 1:length(sampleids)){
        sampleid <- sampleids[k]
        file.footprint <- files.footprint[sampleid]

        # LOAD DATA ---
        dat <- data.table::fread(file=file.footprint, sep="\t", header=TRUE, nThread=1, data.table=FALSE, stringsAsFactors=FALSE, verbose=FALSE)
        colnames(dat)[which(colnames(dat) == "mcrpc_bound")] <- "bound"
        colnames(dat)[which(colnames(dat) == "mcrpc_score")] <- "footprint_score"

        # EXTRACT DATA BY BOUND STATUS ---
        dat <- subset(dat, dat$bound == bound_status)
        
        # PREPARE DATA ---
        dat <- subset(dat, select=c("TFBS_chr","TFBS_start","TFBS_end","TFBS_strand","TFBS_score","peak_chr","peak_start","peak_end","footprint_score"))
        dat$motif_position <- apply(dat, 1, function(x) paste(x[1], as.numeric(x[2]), as.numeric(x[3]), sep="_"))
        dat$peak_position <- apply(dat, 1, function(x) paste(x[6], as.numeric(x[7]), as.numeric(x[8]), sep="_"))
        dat <- subset(dat, select=c("TFBS_chr","TFBS_start","TFBS_end","TFBS_strand","motif_position","peak_position","TFBS_score","footprint_score"))
        colnames(dat) <- c("chr","start","end","strand","motif_position","peak_position","TFBS_score","footprint_score")
        dat$motif_prefix <- motif_prefix
        dat$SampleID <- sampleid

        list.motif_sampleid[[sampleid]] <- dat

        cat(format(Sys.time(), "%b %d %X"), motif_prefix, sampleid, "\n", sep=" ")
    }

    return(list.motif_sampleid)
}


######################################################################
### LOAD ATAC-SEQ MASTER DATA ---
list.rds_atacseq_masterdata <- readRDS(file=file.rds_atacseq_masterdata)

### GET METADATA ---
metadata <- list.rds_atacseq_masterdata$metadata_mcrpc
sample_names <- metadata$Sample_Name

### GET MOTIFS ---
motif_db <- ParseMotifDB(file.motifdb)
motif_db <- getFolderMotifdb(motif_db, dir.diff=file.path(file.path(dir.footprint, sample_names[1]), "bindetect"))
#motif_db <- getFolderMotifdb(motif_db, dir.diff=file.path(file.path(dir.footprint_original, sample_names[1]), "bindetect"))





### EXTRACT MOTIF POSITION ---
list.motif_bound <- list()
for(i in 1:nrow(motif_db)){
    motif_prefix <- motif_db$motif_prefix[i]
    motif_folder <- motif_db$motif_folder[i]

    cat(format(Sys.time(), "%b %d %X"), "START BOUND:", i, motif_prefix, "\n", sep=" ")

        list.motif_bound[[motif_prefix]] <- getMotifPosition(dir.footprint_input=data_footprints, sampleids=sample_names, motif_prefix, motif_folder, bound_status=1)
        #list.motif_bound[[motif_prefix]] <- getMotifPosition(dir.footprint_input=dir.footprint_original, sampleids=sample_names, motif_prefix, motif_folder, bound_status=1)

    cat(format(Sys.time(), "%b %d %X"), "PROCESSED BOUND:", i, motif_prefix, "\n", sep=" ")
}


cat(format(Sys.time(), "%b %d %X"), "WRITING OUTPUT ...", "\n", sep=" ")
    ### SAVE OBJECT TO RDATA FILE ---
    file.rds_motif_bound <- file.path(dir.footprint, "data_footprints/mcrpc_samplewise_motif_position_bound.rds")
    saveRDS(object=list.motif_bound, file=file.rds_motif_bound)
cat(format(Sys.time(), "%b %d %X"), "DONE!", "\n", sep=" ")





### EXTRACT MOTIF POSITION ---
list.motif_unbound <- list()
for(i in 1:nrow(motif_db)){
    motif_prefix <- motif_db$motif_prefix[i]
    motif_folder <- motif_db$motif_folder[i]

    list.motif_unbound[[motif_prefix]] <- getMotifPosition(dir.footprint_input=dir.footprint, sampleids=sample_names, motif_prefix, motif_folder, bound_status=0)
    #list.motif_unbound[[motif_prefix]] <- getMotifPosition(dir.footprint_input=dir.footprint_original, sampleids=sample_names, motif_prefix, motif_folder, bound_status=0)

    cat(format(Sys.time(), "%b %d %X"), "PROCESSED UNBOUND:", i, motif_prefix, "\n", sep=" ")
}

cat(format(Sys.time(), "%b %d %X"), "WRITING OUTPUT ...", "\n", sep=" ")
    ### SAVE OBJECT TO RDATA FILE ---
    file.rds_motif_unbound <- file.path(dir.footprint, "data_footprints/mcrpc_samplewise_motif_position_unbound.rds")
    saveRDS(object=list.motif_unbound, file=file.rds_motif_unbound)
cat(format(Sys.time(), "%b %d %X"), "DONE!", "\n", sep=" ")
