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

dir.pomerantz <- file.path(dir.wrk, "external_data/pomerantz_chipseq_2020")
dir.chipseq_ar <- file.path(dir.wrk, "external_data/chipseq_ar")
dir.chipseq_znf263 <- file.path(dir.wrk, "external_data/znf263_chipseq")


### DEFINE FILES ---
file.pomerantz_mcrpc_AR <- file.path(dir.pomerantz, "pomerantz_chipseq_AR_mcrpc.bed")
file.pomerantz_mcrpc_FOXA1 <- file.path(dir.pomerantz, "pomerantz_chipseq_FOXA1_mcrpc.bed")
file.pomerantz_mcrpc_HOXB13 <- file.path(dir.pomerantz, "pomerantz_chipseq_HOXB13_mcrpc.bed")
file.pomerantz_mcrpc_H3K27ac <- file.path(dir.pomerantz, "pomerantz_chipseq_H3K27ac_mcrpc.bed")

file.chipseq_ar_jin <- file.path(dir.chipseq_ar, "AR_DU145_Chipseq_Jin2014_PMID24875621.bed")
file.chipseq_ar_zhao <- file.path(dir.chipseq_ar, "AR_C42B_Chipseq_Zhao2016_PMID27068475.bed")
file.chipseq_ar_cai <- file.path(dir.chipseq_ar, "AR_22RV1_Chipseq_Cai2018_PMID30270106.bed")

file.chipseq_znf263_imbeault <- file.path(dir.chipseq_znf263, "GSM2466508_Imbeault_ZNF263_chipseq.bed")
file.chipseq_znf263_frietze <- file.path(dir.chipseq_znf263, "GSM476651_Frietze_ZNF263_chipseq.bed")
file.chipseq_znf263_pope <- file.path(dir.chipseq_znf263, "GSM935532_Pope_ZNF263_chipseq.bed")

file.chipseq_myc_barfeld <- file.path(dir.wrk, "external_data/barfeld_chipseq_myc_2017/GSM1907205_Barfeld_MYC_chipseq_peaks_MYC_R1881_DOX_2.bed")
file.chipseq_myc_se <- file.path(dir.wrk, "external_data/se_chipseq_myc_2022/GSE164777_U2OS_MYC_DOX+_peaks_hg38.bed")
file.chipseq_myc_guo <- file.path(dir.wrk, "external_data/guo_chipseq_myc_2021_vcap_GSE157105/GSM4753177_hg38_Peak_VCAP-MYC-0.bed")

file.chipseq_h3k27ac_kron <- file.path(dir.wrk, "external_data/h3k27ac_chipseq_kron_etal_GSE96652/H3K27ac_sum_hg38.bed")

### FUNCTION: loadBED()
loadBED <- function(file.bed){
    # LOAD BED FILES ---
    bed <- data.table::fread(file=file.bed, sep="\t", header=FALSE, nThread=30, data.table=TRUE, verbose=FALSE)
    
    # ADD COLUMN NAMES ---
    if(ncol(bed) == 5){
        colnames(bed) <- c("chr","start","end","peaks","score")
    }else if(ncol(bed) == 4){
        colnames(bed) <- c("chr","start","end","score")        
    }else{
        colnames(bed) <- c("chr","start","end")
    }

    # GRANGES ---
    gr <- GenomicRanges::makeGRangesFromDataFrame(bed, keep.extra.columns=TRUE, ignore.strand=TRUE)

    return(gr)    
}


###################################################################################################################################################
######## POMERANTZ mCRPC ######## 
### COMBINE FILES ---
files.pomerantz_mcrpc <- c(file.pomerantz_mcrpc_AR, file.pomerantz_mcrpc_FOXA1, file.pomerantz_mcrpc_HOXB13, file.pomerantz_mcrpc_H3K27ac)
names(files.pomerantz_mcrpc) <- c("AR","FOXA1","HOXB13","H3K27ac")


### LOAD CHIPSEQ DATA ---
list.pomerantz_mcrpc <- list()
for(i in 1:length(files.pomerantz_mcrpc)){
    protein <- names(files.pomerantz_mcrpc)[i]
    list.pomerantz_mcrpc[[protein]] <- loadBED(file.bed=files.pomerantz_mcrpc[protein])
}

###################################################################################################################################################
######## CHIPSEQ AR ######## 
### COMBINE FILES ---
files.chipseq_ar <- c(file.chipseq_ar_jin, file.chipseq_ar_zhao, file.chipseq_ar_cai)
names(files.chipseq_ar) <- c("Jin2014","Zhao2016","Cai2018")

### LOAD CHIPSEQ DATA ---
list.chipseq_ar <- list()
for(i in 1:length(files.chipseq_ar)){
    dataset <- names(files.chipseq_ar)[i]
    list.chipseq_ar[[dataset]] <- loadBED(file.bed=files.chipseq_ar[dataset])
}



###################################################################################################################################################
######## CHIPSEQ ZNF263 ######## 
### COMBINE FILES ---
files.chipseq_znf263 <- c(file.chipseq_znf263_imbeault, file.chipseq_znf263_frietze, file.chipseq_znf263_pope)
names(files.chipseq_znf263) <- c("Imbeault","Frietze","Pope")


### LOAD CHIPSEQ DATA ---
list.chipseq_znf263 <- list()
for(i in 1:length(files.chipseq_znf263)){
    dataset <- names(files.chipseq_znf263)[i]
    list.chipseq_znf263[[dataset]] <- loadBED(file.bed=files.chipseq_znf263[dataset])
}

###################################################################################################################################################
######## CHIPSEQ MYC ######## 
### COMBINE FILES ---
files.chipseq_myc <- c(file.chipseq_myc_barfeld, file.chipseq_myc_se, file.chipseq_myc_guo)
names(files.chipseq_myc) <- c("Barfeld","Se","Guo")

### LOAD CHIPSEQ DATA ---
list.chipseq_myc <- list()
for(i in 1:length(files.chipseq_myc)){
    dataset <- names(files.chipseq_myc)[i]
    list.chipseq_myc[[dataset]] <- loadBED(file.bed=files.chipseq_myc[dataset])
}



###################################################################################################################################################
######## CHIPSEQ H3K27ac ######## 
gr_chipseq_h3k27ac_kron <- loadBED(file.bed=file.chipseq_h3k27ac_kron)

###################################################################################################################################################
### ADD TO LIST ---
list.output <- list(pomerantz_mcrpc=list.pomerantz_mcrpc,
                    chipseq_ar=list.chipseq_ar,
                    chipseq_znf263=list.chipseq_znf263,
                    chipseq_myc=list.chipseq_myc,
                    chipseq_h3k27ac_kron=gr_chipseq_h3k27ac_kron)


### SAVE OBJECT TO RDATA FILE ---
file.rds_chipseq_bed <- file.path(dir.reproduce_data, "chipseq_bed.rds")
saveRDS(object=list.output, file=file.rds_chipseq_bed)
