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
  optparse::make_option(c("-d", "--dir_input"), type = "character", help = "Path to input file [Required]"),
  optparse::make_option(c("-o", "--file_output"), type = "character", help = "Path to output file [Required]"),
  optparse::make_option(c("-c", "--chr"), type = "character", help = "Chromosome [Required]")
)

parseobj <- optparse::OptionParser(option_list=list.option, usage="usage: Rscript %prog[options]")
opt <- parse_args(parseobj)

### PASS OBJECT VALUES ---
dir.input <- opt$dir_input
file.output <- opt$file_output
chr <- opt$chr

### ITERATE TO PARSE PERM ---
for(k in 1:1000){
    key <- paste(chr, k, sep="_")
    id <- paste("n", k, sep="_")
    file.dat <- file.path(dir.input, paste("perm_", key, ".tsv.gz", sep=""))

    # LOAD DATA ---
    dat <- data.table::fread(file=file.dat, sep="\t", header=TRUE, nThread=5, data.table=FALSE, verbose=FALSE)

    if(k == 1){
        items <- paste("n", c(1:1000), sep="_")
        x <- nrow(dat)
        y <- length(items) + 2
        mat <- matrix(NA, nrow=x, ncol=y, dimnames=list(c(1:x), c("Peak","Gene",items)))
        dt_mat <- data.table::as.data.table( mat )
        mat <- NULL

        # FILL MATRIX ---
        dt_mat$Peak <- dat$Peak
        dt_mat$Gene <- dat$Gene
    }

    # FILL MATRIX ---
    index <- which(colnames(dt_mat) == id)
    dt_mat[,index] <- round(dat$Corr, 4)

    # FLUSH DATA ---
    dat <- NULL
    index <- NULL

    cat(format(Sys.time(), "%b %d %X"), "PROCESSED:", key, "\n", sep=" ")
}


### WRITE OUTPUT ---
cat(format(Sys.time(), "%b %d %X"), "WRITING OUTPUT ... ", "\n", sep=" ")
    write.table(dt_mat, file.output, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
cat(format(Sys.time(), "%b %d %X"), "DONE! ", "\n", sep=" ")


cat(format(Sys.time(), "%b %d %X"), "COMPRESSING ... ", "\n", sep=" ")
    # COMPRESS ---
    cmd <- paste("gzip", file.output, sep=" ")
    system(cmd)    
cat(format(Sys.time(), "%b %d %X"), "DONE! ", "\n", sep=" ")
