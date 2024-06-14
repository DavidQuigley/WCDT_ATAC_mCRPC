###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### LOAD LIBRARIES ---
library("stringr")
library("data.table")

### DEFINE PATH ---
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/analysis/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

dir.tfpeaks2gene <- file.path(dir.reproduce_data, "tfpeaks2gene")
dir.footprints <- file.path(dir.reproduce_data, "tf_footprints/mcrpc_subtypes/data_footprints") 
dir.tf2genes_intermediate <- file.path(dir.reproduce_data, "tf_footprints/mcrpc_subtypes/tf2genes_intermediate") #CREATE


### LOAD TF2GENES FOR EACH MOTIF ---
list.dat <- list()
for(k in 1:458){
    cat(format(Sys.time(), "%b %d %X"), "START:", k, "\n", sep=" ")

    file.dat <- file.path(dir.tf2genes_intermediate, paste("tf2genes_", k, ".tsv.gz", sep=""))

    if( file.exists(file.dat) ){
        list.dat[[k]] <- data.table::fread(file=file.dat, sep="\t", header=TRUE, nThread=50, data.table=FALSE, verbose=FALSE)
    }

    cat(format(Sys.time(), "%b %d %X"), "PROCESSED:", k, "\n", sep=" ")
}


### COMPILE RESULTS ---
df <- do.call(rbind.data.frame, list.dat)
rownames(df) <- NULL


### WRITE OUTPUT ---
cat(format(Sys.time(), "%b %d %X"), "WRITING OUTPUT ... ", "\n", sep=" ")
    file.output <- file.path(dir.tfpeaks2gene, "wcdt_tf2genes_r_0p4_pval_0p05_distTSS_500kb.tsv")
    write.table(df, file.output, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
cat(format(Sys.time(), "%b %d %X"), "DONE! ", file.output, "\n", sep=" ")


### COMPRESS OUTPUT ---
cat(format(Sys.time(), "%b %d %X"), "COMPRESSING FILE ... ", "\n", sep=" ")
    cmd <- paste("gzip", file.output, sep=" ")
    system(cmd)
cat(format(Sys.time(), "%b %d %X"), "DONE! ", "\n", sep=" ")

