###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### LOAD LIBRARIES ---
library("stringr")
library("data.table")
library("dplyr")

### DEFINE PATH ---
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

#dir.footprints <- file.path(dir.wrk, "analysis/12_tf_footprint/subtypes_AR/tobias")
dir.footprints <- file.path(dir.reproduce_data, "tf_footprints/mcrpc_subtypes")
dir.bam <- file.path(dir.reproduce_data, "tf_footprints/mcrpc_subtypes/bam")
dir.bed <- file.path(dir.reproduce_data, "tf_footprints/mcrpc_subtypes/bed")

### DEFINE FILES ---
file.motifdb <- file.path(dir.reproduce, "reference/JASPAR2018_motif_gene_conversion.tsv")
file.rds_footprints <- file.path(dir.footprints, "data_footprints/tobias_mcrpc_subtypes_tf_footprints_overview.rds")

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


### FUNCTION: loadFootprintData() ---
loadFootprintData <- function(file.rds_footprints, grp_id, motif_id){
    # GET FOOTPRINT DATA ---
    dat <- list.rds_footprints[[grp_id]][[motif_id]]
    dat$FeatureID <- apply(dat, 1, function(x) paste(x[1],as.numeric(x[2]),as.numeric(x[3]), sep="_"))

    # TRIM DATA BY COLUMNS ---
    dat <- subset(dat, select=c("FeatureID",colnames(dat)[stringr::str_detect(colnames(dat), "_bound")]) )

    return(dat)
}


### FUNCTION: getFootprintsByState() ---
getFootprintsByState <- function(dat, grp_id, state){
    # GROUP IDS ---
    grp_id_1 <- stringr::str_split(grp_id, "_")[[1]][1]
    grp_id_2 <- stringr::str_split(grp_id, "_")[[1]][2]

    # COLUMN INDEX ---
    index_1 <- which( colnames(dat) == paste(grp_id_1, "bound", sep="_") )
    index_2 <- which( colnames(dat) == paste(grp_id_2, "bound", sep="_") )

    # GET FEATURES ----
    featureid_1 <- dat$FeatureID[which(dat[,index_1] == state)]
    featureid_2 <- dat$FeatureID[which(dat[,index_2] == state)]

    # STORE DATA ---
    list.features <- list()
    list.features[[grp_id_1]] <- featureid_1
    list.features[[grp_id_2]] <- featureid_2

    return(list.features)
}


### FUNCTION: getFootprintsAggregate() ---
getFootprintsAggregate <- function(list.features, grp_ids, subtypes){
    list.features_subtype <- list()
    for(subtype in subtypes){
        g_ids <- grp_ids[stringr::str_detect(grp_ids, subtype)]

        list_g <- list()
        for(itr in 1:length(g_ids)){
            y <- unlist(list.features[[ g_ids[itr] ]][ subtype ])
            names(y) <- NULL 
            list_g[[itr]] <- y
        }

        list.features_subtype[[subtype]] <- unique(unlist(list_g))
    }

    return(list.features_subtype)
}



#########################################################################################################################
### LOAD FOOTPRINTS ---
list.rds_footprints <- readRDS(file=file.rds_footprints)

### GET MOTIFS ---
motif_db <- ParseMotifDB(file.motifdb)
motif_db <- subset(motif_db, motif_db$motif_id %in% names( list.rds_footprints[[1]] ) )
rownames(motif_db) <- NULL


#########################################################################################################################
### COMPILE DATA FOR EACH MOTIF --------------
list.footprints_bound <- list()
list.footprints_unbound <- list()
n <- nrow(motif_db)
for(i in 1:n){
    motif_prefix <- motif_db$motif_prefix[i]
    motif_name <- motif_db$motif_name[i]
    motif_id <- motif_db$motif_id[i]

    cat(format(Sys.time(), "%X"), "START:", i, "OF", n, motif_prefix, "\n", sep=" ")  

    grp_ids <- dm_grp$id
    keys <- paste(grp_ids, motif_prefix, sep="_")

    # GET TOBIAS SUMMARY DATA ---
    list.features_bound <- list()
    list.features_unbound <- list()
    for(j in 1:length(keys)){
        grp_id <- grp_ids[j]
        key <- keys[j]

        # GET TF BINDING SITES PER GROUP ---
        dat <- loadFootprintData(file.rds_footprints, grp_id, motif_id)
        list.features_bound[[grp_id]] <- getFootprintsByState(dat, grp_id, state=1)
        list.features_unbound[[grp_id]] <- getFootprintsByState(dat, grp_id, state=0)
    }

    # GET FOOTPRINTS AGGREGATE ----
    list.footprints_bound[[motif_prefix]] <- getFootprintsAggregate(list.features=list.features_bound, grp_ids, subtypes)
    list.footprints_unbound[[motif_prefix]] <- getFootprintsAggregate(list.features=list.features_unbound, grp_ids, subtypes)

    cat(format(Sys.time(), "%X"), "DONE:", i, "OF", n, motif_prefix, "\n", sep=" ")   
}




###################################################################################################################################################
### ADD TO LIST ---
list.output <- list(bound=list.footprints_bound,
                    unbound=list.footprints_unbound)


### SAVE OBJECT TO RDATA FILE ---
file.rds_tf_footprints_bound_unbound <- file.path(dir.footprints, "data_footprints/tobias_mcrpc_subtypes_tf_footprints_featureids_bound_unbound.rds")
saveRDS(object=list.output, file=file.rds_tf_footprints_bound_unbound)
