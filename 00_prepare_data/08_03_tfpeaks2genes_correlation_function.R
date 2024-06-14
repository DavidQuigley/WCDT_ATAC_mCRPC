###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### LOAD LIBRARIES ---
suppressPackageStartupMessages(require("stringr", quietly=TRUE, warn.conflicts=FALSE))
suppressPackageStartupMessages(require("data.table", quietly=TRUE, warn.conflicts=FALSE))
suppressPackageStartupMessages(require("optparse", quietly=TRUE, warn.conflicts=FALSE))

### PASS ARGUMENTS ---
list.option <- list(
  optparse::make_option(c("-a", "--file_atacseq"), type = "character", help = "Path to input atacseq file [Required]"),
  optparse::make_option(c("-r", "--file_rnaseq"), type = "character", help = "Path to input rnaseq file [Required]"),
  optparse::make_option(c("-o", "--file_output"), type = "character", help = "Path to output file [Required]")
)

parseobj <- optparse::OptionParser(option_list=list.option, usage="usage: Rscript %prog[options]")
opt <- parse_args(parseobj)

### PASS OBJECT VALUES ---
file.atacseq <- opt$file_atacseq
file.rnaseq <- opt$file_rnaseq
file.output <- opt$file_output


### FUNCTION: matcor() ---
# FROM CCA R-package
matcor <- function (X, Y) {
    matcorXY = cor(cbind(X, Y), use="pairwise", method="spearman")
    return(list(XYcor = matcorXY))
}

### FUNCTION: loadData() ---
loadData <- function(file.dat){
    dat <- data.table::fread(file=file.dat, sep="\t", header=TRUE, nThread=10, data.table=FALSE, verbose=FALSE)
    rownames(dat) <- dat$V1
    dat$V1 <- NULL
    mat <- as.matrix(dat)
    return(mat)
}


### LOAD DATA ---
cat(format(Sys.time(), "%b %d %X"), "LOAD DATA ... ", "\n", sep=" ")
    mat_atacseq <- loadData(file.dat=file.atacseq)
    mat_rnaseq <- loadData(file.dat=file.rnaseq)
cat(format(Sys.time(), "%b %d %X"), "DONE! ", "\n", sep=" ")

### GET CORRELATION ---
cat(format(Sys.time(), "%b %d %X"), "START: CORRELATION ... ", "\n", sep=" ")
    list.res_cor <- matcor(X=mat_atacseq, Y=mat_rnaseq)
cat(format(Sys.time(), "%b %d %X"), "COMPLETE! ", "\n", sep=" ")

### GET CORRELATION MATRIX ---
dat <- list.res_cor$XYcor

### PARSE PAIR-WISE CORRELATION MATRIX ---
cat(format(Sys.time(), "%b %d %X"), "PARSE CORRELATION MATRIX ...", "\n", sep=" ")
# GET INDICES ---
index.atacseq <- which( stringr::str_detect(rownames(dat), "_") == TRUE )
rnaseq_start <- max(index.atacseq) + 1
rnaseq_end <- length( rownames(dat) )
index.rnaseq <- c(rnaseq_start:rnaseq_end)

### EXTRACT MATRIX ---
mat <- as.matrix(dat)
rownames(mat) <- rownames(dat)
mat <- mat[index.atacseq, index.rnaseq]

dat <- NULL

### RESHAPE MATRIX ---
dm <- reshape2::melt(mat)
colnames(dm) <- c("Peak","Gene","Corr")
dm$Corr <- round(dm$Corr, 4)
    
cat(format(Sys.time(), "%b %d %X"), "DONE! ", "\n", sep=" ")


### WRITE OUTPUT ---
cat(format(Sys.time(), "%b %d %X"), "WRITING OUTPUT ... ", "\n", sep=" ")
    write.table(dm, file.output, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
cat(format(Sys.time(), "%b %d %X"), "DONE! ", "\n", sep=" ")

cat(format(Sys.time(), "%b %d %X"), "COMPRESSING ... ", "\n", sep=" ")
    # COMPRESS ---
    cmd <- paste("gzip", file.output, sep=" ")
    system(cmd)    
cat(format(Sys.time(), "%b %d %X"), "DONE! ", "\n", sep=" ")
