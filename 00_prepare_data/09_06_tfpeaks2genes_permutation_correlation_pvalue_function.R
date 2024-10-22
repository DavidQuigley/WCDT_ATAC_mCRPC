###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### LOAD LIBRARIES ---
suppressPackageStartupMessages(require("stringr", quietly=TRUE, warn.conflicts=FALSE))
suppressPackageStartupMessages(require("data.table", quietly=TRUE, warn.conflicts=FALSE))
suppressPackageStartupMessages(require("foreach", quietly=TRUE, warn.conflicts=FALSE))
suppressPackageStartupMessages(require("doParallel", quietly=TRUE, warn.conflicts=FALSE))
suppressPackageStartupMessages(require("zoo", quietly=TRUE, warn.conflicts=FALSE))
suppressPackageStartupMessages(require("optparse", quietly=TRUE, warn.conflicts=FALSE))


### PASS ARGUMENTS ---
list.option <- list(
  optparse::make_option(c("-c", "--file_corr"), type = "character", help = "Path to input correlation file [Required]"),
  optparse::make_option(c("-p", "--file_permute"), type = "character", help = "Path to input permutation file [Required]"),
  optparse::make_option(c("-o", "--file_output"), type = "character", help = "Path to output file [Required]")
)

parseobj <- optparse::OptionParser(option_list=list.option, usage="usage: Rscript %prog[options]")
opt <- parse_args(parseobj)

### PASS OBJECT VALUES ---
file.corr <- opt$file_corr
file.permute <- opt$file_permute
file.output <- opt$file_output


###############################################################################################################
#### FUNCTION: densityauc() ------
#' EXTRACTED FROM vulcan R package https://www.bioconductor.org/packages/release/bioc/html/vulcan.html
#' densityauc - Calculate the AUC of a density object
#'
#' This function will calculate the AUC of a density object generated by the
#' \code{'density'} function.
#'
#' @param dens a density object
#' @param window a vector with two values, specifying the left and right borders
#' for the AUC to be calculated
#' @return a numeric value for the density AUC
#' @examples
#' set.seed(1)
#' a<-rnorm(1000)
#' d<-density(a)
#' window<-c(2,3)
#' da<-densityauc(d,window)
#'
#' plot(d,main='')
#' abline(v=window,lty=2)
#' title(paste0('AUC between lines=',da))
#'
#'
#'
#' @export
densityauc <- function(dens, window){
    xt <- diff(dens$x[dens$x > window[1] & dens$x < window[2]])
    yt <- zoo::rollmean(dens$y[dens$x > window[1] & dens$x < window[2]], 2)
    return( sum(xt * yt) )
}
###############################################################################################################


### FUNCTION: compute_pvalue() ---
compute_pvalue <- function(val_permute, val_obs){
    # COMPUTE DENSITY OF RANDOM (PERMUTED) DISTRIBUTION ----
    dens <- density(val_permute)

    # COMPUTE AUC ---
    if(val_obs <= mean(val_permute)){
        y_auc <- densityauc(dens, window=c(min(val_permute), val_obs))
    }else{
        y_auc <- densityauc(dens, window=c(val_obs, max(val_permute)))
    }

    # COMPUTE PVALUE (SUM OF LEFT AND RIGHT AUC) ----
    pval <- y_auc * 2

    return(pval)
}


### FUNCTION: getPvalues() ---
getPvalues <- function(dat_corr, mat){
    n <- nrow(dat_corr)
    for(k in 1:n){
        cat(format(Sys.time(), "%b %d %X"), "START:", k, "OF", n, "\n", sep=" ")

        val_obs <- dat_corr$Corr[k]
        val_permute <- as.numeric(mat[k,])

        # COMPUTE PVALUE ---
        dat_corr$pvalue[k] <- compute_pvalue(val_permute, val_obs)
        
        cat(format(Sys.time(), "%b %d %X"), "PROCESSED:", k, "OF", n, "\n", sep=" ")
    }

    return(dat_corr)
}

### FUNCTION: getPvalues() ---
getPvalues_parallel <- function(dat_corr, mat){
    n <- nrow(dat_corr)

    # DECLARE CLUSTER ----
	no_cores <- 55
	cl <- makeCluster(no_cores)
	registerDoParallel(cl)

    lpar <- foreach(k = 1:n, .inorder=FALSE, .combine='c', .errorhandling="pass", .export=c("compute_pvalue","densityauc"), .packages="zoo", .verbose=TRUE) %dopar%{	
    #for(k in 1:n){
        #cat(format(Sys.time(), "%b %d %X"), "START:", k, "OF", n, "\n", sep=" ")

        val_obs <- dat_corr$Corr[k]
        val_permute <- as.numeric(mat[k,])

        # COMPUTE PVALUE ---
        #dat_corr$pvalue[k] <- compute_pvalue(val_permute, val_obs)
        lpar <- compute_pvalue(val_permute, val_obs)
        
        #cat(format(Sys.time(), "%b %d %X"), "PROCESSED:", k, "OF", n, "\n", sep=" ")
    }

    stopCluster(cl)

    # ADD PVALUES ---
    dat_corr$pvalue <- lpar

    return(dat_corr)
}


### FUNCTION: partitionData() ---
partitionData <- function(dat){
    n <- 1000
    n1 <- n + 1
    q <- nrow(dat) %/% n
    x <- 1

    list.dat <- list()
    for(k in 1:n1){
        cat(format(Sys.time(), "%b %d %X"), "PARTITION DATA START:", k, "OF", n1, "\n", sep=" ")

        # GET INDEX RANGE ----
        y <- x + (q-1)

        if(y > nrow(dat)){
            y <- nrow(dat)
        }

        # EXTRACT DATA ---
        list.dat[[k]] <- dat[x:y,]

        x <- y + 1

        cat(format(Sys.time(), "%b %d %X"), "PARTITION DATA STOP:", k, "OF", n1, "\n", sep=" ")
    }

    return(list.dat)
}

### FUNCTION: getPvaluesPartition() ---
getPvaluesPartition <- function(list.corr, list.perm){
    list.pval <- list()
    n <- 1000
    n1 <- n + 1
    for(k in 1:n1){
        cat(format(Sys.time(), "%b %d %X"), "START:", k, "\n", sep=" ")

        list.pval[[k]] <- getPvalues_parallel(dat_corr=list.corr[[k]], mat=list.perm[[k]])

        cat(format(Sys.time(), "%b %d %X"), "PROCESSED ", k, "\n", sep=" ")
    }
    return(list.pval)
}


### FUNCTION: getPermutationPvalues() ----
getPermutationPvalues <- function(file.corr, file.permute){
    cat(format(Sys.time(), "%b %d %X"), "LOADING INPUT DATA ... ", "\n", sep=" ")
        # LOAD DATA ---
        dat_permute <- data.table::fread(file=file.permute, sep="\t", header=TRUE, nThread=20, data.table=FALSE, verbose=FALSE)
        dat_corr <- data.table::fread(file=file.corr, sep="\t", header=TRUE, nThread=10, data.table=FALSE, verbose=FALSE)
        dat_corr$pvalue <- NA
    cat(format(Sys.time(), "%b %d %X"), "DONE! ", "\n", sep=" ")

    # MATCH DATA ---
    y_perm <- apply(dat_permute, 1, function(x) paste(x[1], x[2], sep=":"))
    y_corr <- apply(dat_corr, 1, function(x) paste(x[1], x[2], sep=":"))
    dat_permute <- dat_permute[match(y_corr, y_perm),]

    y_perm <- NULL
    y_corr <- NULL

    # GET MATRIX ---
    mat <- as.matrix(dat_permute[,3:ncol(dat_permute)])
    dat_permute <- NULL

    # PARTITION DATA ---
    list.corr <- partitionData(dat=dat_corr)
    list.perm <- partitionData(dat=mat)
    
    # FLUSH DATA ---
    dat_corr <- NULL
    mat <- NULL

    # GET PVALUES ---
    list.pval <- getPvaluesPartition(list.corr, list.perm)
    
    # FLUSH DATA ---
    list.corr <- NULL
    list.perm <- NULL

    # AGGREGATE DATA ---
    df <- do.call(rbind.data.frame, list.pval)
    rownames(df) <- NULL

    return(df)
}



### GET PERMUTATION PVALUES ---
df <- getPermutationPvalues(file.corr, file.permute)

### WRITE OUTPUT ---
cat(format(Sys.time(), "%b %d %X"), "WRITING OUTPUT ... ", "\n", sep=" ")
    write.table(df, file.output, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
cat(format(Sys.time(), "%b %d %X"), "DONE! ", "\n", sep=" ")

### COMPRESS OUTPUT ---
cat(format(Sys.time(), "%b %d %X"), "COMPRESSING FILE ... ", "\n", sep=" ")
    cmd <- paste("gzip", file.output, sep=" ")
    system(cmd)
cat(format(Sys.time(), "%b %d %X"), "DONE! ", "\n", sep=" ")

