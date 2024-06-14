###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### LOAD LIBRARIES ---
suppressPackageStartupMessages(require("stringr", quietly=TRUE, warn.conflicts=FALSE))
suppressPackageStartupMessages(require("data.table", quietly=TRUE, warn.conflicts=FALSE))
suppressPackageStartupMessages(require("optparse", quietly=TRUE, warn.conflicts=FALSE))

#######################################################################################################################################
### PASS ARGUMENTS ---
list.option <- list(
    optparse::make_option(c("-t", "--file_tf2peaks"), type = "character", help = "Path to input tf2peaks file [Required]"),
    optparse::make_option(c("-p", "--file_peaks2gene"), type = "character", help = "Path to input peaks2gene file [Required]"),
    optparse::make_option(c("-m", "--motif_id"), type = "character", help = "motif [Required]"),
    optparse::make_option(c("-o", "--file_output"), type = "character", help = "Path to output file [Required]")
)

parseobj <- optparse::OptionParser(option_list=list.option, usage="usage: Rscript %prog[options]")
opt <- parse_args(parseobj)

### PASS OBJECT VALUES ---
file.tf2peaks <- opt$file_tf2peaks
file.peaks2gene <- opt$file_peaks2gene
motif_prefix <- opt$motif_id
file.output <- opt$file_output


#######################################################################################################################################
### GET TF2PEAKS DATA ---
list.motif_footprints <- readRDS(file=file.tf2peaks)
list.motif_footprints_subtypes <- list.motif_footprints[[motif_prefix]]
dat_tf2peaks <- do.call(rbind.data.frame, list.motif_footprints_subtypes)
rownames(dat_tf2peaks) <- NULL
dat_tf2peaks <- subset(dat_tf2peaks, select=c("motif_prefix","peak_position"))
n <- nrow(dat_tf2peaks)



#######################################################################################################################################
### GET PEAKS2GENE DATA ---
dat_peaks2gene <- data.table::fread(file=file.peaks2gene, sep="\t", header=TRUE, nThread=10, data.table=FALSE, verbose=FALSE)
dat_peaks2gene <- subset(dat_peaks2gene, select=c("Peak","Gene","Corr"))

### ASSOCIATE TF2GENES ---
list.df <- list()    
for(i in 1:n){
    val_peak <- dat_tf2peaks$peak_position[i]
    index <- which(dat_peaks2gene$Peak == val_peak)

    if(length(index) != 0){
        dat_temp <- dat_peaks2gene[index,]
        list.df[[i]] <- cbind(motif_prefix=motif_prefix, dat_temp)
    }

    cat(format(Sys.time(), "%b %d %X"), "PROCESSED:", i, "OF", n, "\n", sep=" ")
}

df <- do.call(rbind.data.frame, list.df)
rownames(df) <- NULL



#######################################################################################################################################
### WRITE OUTPUT ---
if(nrow(df) != 0){
    cat(format(Sys.time(), "%b %d %X"), "WRITING OUTPUT ... ", "\n", sep=" ")
        write.table(df, file.output, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
    cat(format(Sys.time(), "%b %d %X"), "DONE! ", "\n", sep=" ")


    # COMPRESS OUTPUT ---
    cat(format(Sys.time(), "%b %d %X"), "COMPRESSING FILE ... ", "\n", sep=" ")
        cmd <- paste("gzip", file.output, sep=" ")
        system(cmd)
    cat(format(Sys.time(), "%b %d %X"), "DONE! ", "\n", sep=" ")
}else{
    cat(format(Sys.time(), "%b %d %X"), "ALERT! NO MATCH FOUND!", "\n", sep=" ")
}

