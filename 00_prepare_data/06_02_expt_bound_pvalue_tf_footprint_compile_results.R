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



### SUBTYPES ---
subtypes <- c("ARpNEn","ARpNEp","ARlNEn","ARnNEp","ARnNEn")

### PVALUES ---
pvals <- c(0.5, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001)
names(pvals) <- c("p5","p1","p05","p01","p005","p001","p0005","p0001")




### GET PAIRWISE ALL GROUP COMBINATIONS ---
dm_grp <- do.call(rbind.data.frame, combn(x=subtypes, m=2, simplify=FALSE)) %>%
                dplyr::mutate_all(as.character)
colnames(dm_grp) <- c("Group1","Group2")
rownames(dm_grp) <- NULL
dm_grp$Analysis <- paste(dm_grp$Group1, dm_grp$Group2, sep="_")
dm_grp <- cbind(dm_grp, cbind(sapply(names(pvals), function(x) { paste(sprintf("%s",x), dm_grp$Analysis, sep="_") } )))

### RESHAPE DATA ---
des <- reshape2::melt(dm_grp, id.vars=c("Group1","Group2","Analysis"), variable.name="pval", value.name="fname")
des$pval <- as.character(des$pval)
des$dir_analysis <- file.path(dir.footprints, des$fname)



### FUNCTION: getTFdata() ---
getTFdata <- function(motifs, dir_analysis){
    files.dat <- sapply(motifs, function(x) { paste(dir_analysis,"/", sprintf("%s/%s",x,x), "_overview.txt", sep="") } )

    # LOOP FOR EACH MOTIF ---
    list.dat <- list()
    for(k in 1:length(motifs)){
        motif <- motifs[k]
        file.dat <- files.dat[motif]
        dat <- data.table::fread(file.dat, data.table=FALSE, header=TRUE, stringsAsFactors=FALSE, nThread=1)   
        list.dat[[motif]] <- dat
    }

    return(list.dat)
}




### LOOP FOR EACH GROUP ---
list.df <- list()
for(i in 1:length(pvals)){
    pval <- names(pvals)[i]

    # SUBSET BY PVALUE GROUP ---
    des_temp <- des[which(des$pval == pval),]
    
    # LOOP FOR EACH ANALYSIS ---
    list.dat_analysis <- list()
    for(j in 1:nrow(des_temp)){
        analysis <- des_temp$Analysis[j]
        dir_analysis <- des_temp$dir_analysis[j]
        list.dat_analysis[[analysis]]<- getTFdata(motifs, dir_analysis)
    }
    
    # STORE DATA ---
    list.df[[pval]] <- list.dat_analysis

    cat("PROCESSED:", pval, "\n", sep=" ")
}


### SAVE OBJECT TO RDATA FILE ---
file.rds_expt_bound_pvalue_footprint_summary <- file.path(dir.footprints, "data_footprints/expt_bound_pvalue_footprint_summary.rds")
saveRDS(object=list.df, file=file.rds_expt_bound_pvalue_footprint_summary)
