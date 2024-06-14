###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### LOAD LIBRARIES ---
cat(format(Sys.time(), "%a %b %d %X %Y"), "LOADING PACKAGE: stringr", "\n", sep=" ")
suppressPackageStartupMessages(require("stringr", quietly=TRUE, warn.conflicts=FALSE))

cat(format(Sys.time(), "%a %b %d %X %Y"), "LOADING PACKAGE: data.table", "\n", sep=" ")
suppressPackageStartupMessages(require("data.table", quietly=TRUE, warn.conflicts=FALSE))

cat(format(Sys.time(), "%a %b %d %X %Y"), "LOADING PACKAGE: GenomicRanges", "\n", sep=" ")
suppressPackageStartupMessages(require("GenomicRanges", quietly=TRUE, warn.conflicts=FALSE))

cat(format(Sys.time(), "%a %b %d %X %Y"), "LOADING PACKAGE: IRanges", "\n", sep=" ")
suppressPackageStartupMessages(require("IRanges", quietly=TRUE, warn.conflicts=FALSE))

cat(format(Sys.time(), "%a %b %d %X %Y"), "LOADING PACKAGE: plyranges", "\n", sep=" ")
suppressPackageStartupMessages(require("plyranges", quietly=TRUE, warn.conflicts=FALSE))

cat(format(Sys.time(), "%a %b %d %X %Y"), "LOADING PACKAGE: ChIPseeker", "\n", sep=" ")
suppressPackageStartupMessages(require("ChIPseeker", quietly=TRUE, warn.conflicts=FALSE))

cat(format(Sys.time(), "%a %b %d %X %Y"), "LOADING PACKAGE: ChIPQC", "\n", sep=" ")
suppressPackageStartupMessages(require("ChIPQC", quietly=TRUE, warn.conflicts=FALSE))

cat(format(Sys.time(), "%a %b %d %X %Y"), "LOADING PACKAGE: TxDb.Hsapiens.UCSC.hg38.knownGene", "\n", sep=" ")
suppressPackageStartupMessages(require("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly=TRUE, warn.conflicts=FALSE))

cat(format(Sys.time(), "%a %b %d %X %Y"), "LOADING PACKAGE: BSgenome.Hsapiens.UCSC.hg38", "\n", sep=" ")
suppressPackageStartupMessages(require("BSgenome.Hsapiens.UCSC.hg38", quietly=TRUE, warn.conflicts=FALSE))

cat(format(Sys.time(), "%a %b %d %X %Y"), "LOADING PACKAGE: aroma.light", "\n", sep=" ")
suppressPackageStartupMessages(require("aroma.light", quietly=TRUE, warn.conflicts=FALSE))


##################################################################################
cat(format(Sys.time(), "%a %b %d %X %Y"), "LOADING FUNCTIONS ...", "\n", sep=" ")
##################################################################################


### FUNCTION: parseGTF() --
#' parseGTF(): Parse GTF files parseGTF() takes in GTF annotation file and parses gene/transcript/exon data
#' 
#' Function to parse GTF annotation file and returns gene/transcript/exon data
#'
#' @param file.gtf path to the GTF file
#' @param feature.type type of feature: "gene", "exon", or "transcript"
#'
#' @return a dataframe containing the feature annotation
#'
#' @author Raunak Shrestha
#'
#' @export
#' FUNCTION: parseGTF() --
parseGTF <- function(file.gtf, feature.type){
    # CHROMOSOMES ---
    chromosome <- paste("chr", c(1:22,"X","Y"), sep="")

    # LOAD GTF FILE ---
    #annot <- data.table::fread(file.gtf, data.table=FALSE, stringsAsFactors=FALSE, skip=5)
    annot <- data.table::fread(file=file.gtf, sep="\t", header=FALSE, nThread=50, data.table=TRUE, stringsAsFactors=FALSE, skip=5, verbose=FALSE)
    annot <- subset(annot, annot$V1 %in% chromosome)

    # PARSE GTF COLUMNS ---
    if(feature.type == "gene"){    
        annot <- subset(annot, annot$V3 == "gene")
    }else if(feature.type == "exon"){
        annot <- subset(annot, annot$V3 == "exon")
    }else if(feature.type == "transcript"){
        annot <- subset(annot, annot$V3 == "transcript")
    }

    list.annot <- noquote(stringr::str_split(annot$V9, "; "))
    annot$EnsemblID <- unlist(lapply(stringr::str_split(unlist(lapply(list.annot, function(x) gsub('"', '', x[stringr::str_detect(x, "gene_id")]))), " "), function(x) x[2]))
    annot$Gene <- unlist(lapply(stringr::str_split(unlist(lapply(list.annot, function(x) gsub('"', '', x[stringr::str_detect(x, "gene_name")]))), " "), function(x) x[2]))
    annot$GeneType <- unlist(lapply(stringr::str_split(unlist(lapply(list.annot, function(x) gsub('"', '', x[stringr::str_detect(x, "gene_type")]))), " "), function(x) x[2]))

    # TRIM DATA ---
    if(feature.type == "gene"){ 
        annot <- subset(annot, annot$GeneType == "protein_coding")   
        items <- c("V1","V4","V5","Gene","EnsemblID","V7")
        df <- subset(annot, select=items)
        colnames(df) <- c("Chrom","Start","End","Gene","EnsemblID","Strand")
    }else if(feature.type == "exon"){
        annot <- subset(annot, annot$GeneType == "protein_coding")
        list.annot <- stringr::str_split(annot$V9, "; ")
        annot$TranscriptID <- unlist(lapply(stringr::str_split(unlist(lapply(list.annot, function(x) gsub('"', '', x[stringr::str_detect(x, "transcript_id")]))), " "), function(x) x[2]))
        annot$ExonID <- unlist(lapply(stringr::str_split(unlist(lapply(list.annot, function(x) gsub('"', '', x[stringr::str_detect(x, "exon_id")]))), " "), function(x) x[2]))  
        items <- c("V1","V4","V5","Gene","EnsemblID","TranscriptID","ExonID","V7")
        df <- subset(annot, select=items)
        colnames(df) <- c("Chrom","Start","End","Gene","GeneID","TranscriptID","ExonID","Strand")
    }else if(feature.type == "transcript"){
        annot$TranscriptID <- unlist(lapply(stringr::str_split(unlist(lapply(list.annot, function(x) gsub('"', '', x[stringr::str_detect(x, "transcript_id")]))), " "), function(x) x[2]))
        items <- c("V1","V4","V5","Gene","EnsemblID","V7","TranscriptID","GeneType")
        df <- subset(annot, select=items)
        colnames(df) <- c("Chrom","Start","End","Gene","EnsemblGeneID","Strand","EnsemblTranscriptID","GeneType")
    }        

    return(df)
}


### FUNCTION: getPeakFiles() ---
getPeakFiles <- function(dir.process, peak_type, sampleid, peaks_shift=TRUE){
    if(peaks_shift == TRUE){
        dir.peak <- file.path( file.path(dir.process, sampleid), "peaks")
    }else {
        dir.peak <- file.path( file.path(dir.process, sampleid), "peaks_uncorrected")
    }
    
    filename <- paste(sampleid, ".filtered.", peak_type, "Peak", sep="")
    file.peak <- file.path(dir.peak, filename)
    
    return(file.peak)
}

### FUNCTION: readPeakFile() ---
readPeakFile <- function(file.peak, peak_type){
    # CHROMOSOME ---
    chromosomes <- paste("chr", c(1:22,"X","Y"), sep="")

    # LOAD PEAK FILE ---
    peak.bed <- read.delim(file.peak, header=FALSE, stringsAsFactors=FALSE)

    # DROP THE 6TH COLUMN ---
    peak.bed$V6 <- NULL

    # ADD COLUMN NAMES ---
    colnames(peak.bed) <- c("chromosome", "start", "end", "peak_id", "int_qvalue", "fold_enrichment_summit", "pvalue", "qvalue")

    if(peak_type == "narrow") {
      colnames(peak.bed)[9] <- "summit_position_to_peak_start"
    }   

    # SUBSET BY CHROMOSOME ---
    peak.bed <- peak.bed[which(peak.bed$chromosome %in% chromosomes),]

    return(peak.bed)
}

### FUNCTION: getPeakBED() ---
getPeakBED <- function(sampleids, dir.process, peak_type, peaks_shift=TRUE){
    list.peaks <- list()

    # GET PEAK PER SAMPLEID ---
    for(sampleid in sampleids){
        # LOAD PEAK FILE ---
        file.peak <- getPeakFiles(dir.process, peak_type, sampleid, peaks_shift)

        if(length(readLines(paste0(file.peak))) != 0){
            dat.peak <- readPeakFile(file.peak, peak_type)
            
            # ADD PEAK ID ---
            dat.peak$peak_id <- paste(sampleid, dat.peak$peak_id, sep="_")
            dat.peak$SampleID <- sampleid
        }else {
            dat.peak <- data.frame()
        }

        # STORE BED ---
        list.peaks[[sampleid]] <- dat.peak

        cat(format(Sys.time(), "%a %b %d %X %Y"), "PROCESSED:", sampleid, "\n", sep=" ")
    }

    return(list.peaks)
}

### FUNCTION: getPeakBEDList() ---
getPeakBEDList <- function(sampleids, dir.process, peak_type, peaks_shift){
    # CHROMOSOME ---
    chromosomes <- paste("chr", c(1:22,"X","Y"), sep="")

    # GET PEAKS ---
    list.peaks <- getPeakBED(sampleids, dir.process, peak_type, peaks_shift)

    # REMOVE DUPLICATE PEAKS ---
    list.atacseq.peaks <- list()

    for(sampleid in sampleids){
        cat("START:", sampleid, "\n", sep="\t")

        # GET PEAKS ---
        dat <- list.peaks[[sampleid]]
        dat <- subset(dat, select=c("chromosome","start","end"))
        dat <- dat[!duplicated(dat),]
        dat <- dat[which(dat$chromosome %in% chromosomes),]
        dat$SampleID <- sampleid
        #dat$FeatureID <- apply(dat, 1, function(x) paste(x[4], x[1], as.numeric(x[2]), as.numeric(x[3]), sep="_"))

        # CONVERT TO GENOMEIC RANGES ---
        gr <- GenomicRanges::makeGRangesFromDataFrame(df=dat, keep.extra.columns=TRUE)
        gr <- sort(gr)

        # CONVERT BACK TO DATAFRAME ---
        df <- as.data.frame(gr)
        df <- subset(df, select=c("seqnames","start","end","SampleID"))
        colnames(df) <- c("chr","start","end","SampleID")
        rownames(df) <- NULL

        # STORE DATA TO LIST ---
        list.atacseq.peaks[[sampleid]] <- df

        cat("DONE:", sampleid, "\n", sep="\t")
    }

    return(list.atacseq.peaks)
}


### FUNCTION: getNonOverlapPeaks() ---
# Compile all ATAC-seq peaks from all samples and then find a unique list of non-overlapping genomic regions
getNonOverlapPeaks <- function(list.peaks, sampleids){
    list.peaks_gr <- lapply(list.peaks, ChIPQC:::GetGRanges, simple = TRUE)
    names(list.peaks_gr) <- sampleids
    
    list.peaks_gr <- GenomicRanges::GRangesList(list.peaks_gr)   
    reduced <- GenomicRanges::reduce(unlist(list.peaks_gr))
    consensusIDs <- paste("consensus", seq(1, length(reduced)), sep="_")
    mcols(reduced) <- do.call(cbind, lapply(list.peaks_gr, function(x) (reduced %over% x) + 0))

    reducedConsensus <- reduced
    mcols(reducedConsensus) <- cbind(as.data.frame(mcols(reducedConsensus)), consensusIDs)
    consensusIDs <- paste("consensus", seq(1, length(reducedConsensus)), sep="_")

    return(reducedConsensus)
}



### FUNCTION: binGenome() ---
binGenome <- function(genome, df.encode.blacklist=NULL, chromosomes, tile.size){
    # ENCODE BLACKLIST
    #gr_encode.blacklist <- GenomicRanges::makeGRangesFromDataFrame(df.encode.blacklist, keep.extra.columns=FALSE)

    # CREATE BIN ---
    bins <- GenomicRanges::tileGenome(seqinfo(genome), tilewidth=tile.size, cut.last.tile.in.chrom=TRUE)

    # FILTER OUT ENCODE BLACKLISTED REGIONS ---
    #bins <- setdiff(bins, gr_encode.blacklist)

    # RESTRICT TO STANDARD CHROMOSOMES ---
    bins <- bins[seqnames(bins) %in% chromosomes]

    return(bins)
}



### FUNCTION: getPeakRegionsFeatureCount() ---
getPeakRegionsFeatureCount <- function(dir.process, list.peaks, sampleids, nthreads){
    # GET NON-OVARLAPPING PEAKS ---
    consensusToCount <- getNonOverlapPeaks(list.peaks, sampleids)


    # RETAIN PEAKS THAT ARE DETEDCTED IN AT LEAST 2 SAMPLES ---
    occurrences <- GenomicRanges::elementMetadata(consensusToCount) %>% 
                                        as.data.frame %>% 
                                        dplyr::select(-consensusIDs) %>% 
                                        rowSums

    consensusToCount <- consensusToCount[occurrences >= 2,]

    # GET ATAC-SEQ PEAK FEATURE COUNT ---
    dat.count <- getATACseqFeatureCount(dir.process, consensusToCount, sampleids, nthreads)

    return(dat.count)
}


### FUNCTION: getATACseqFeatureCount() ---
getATACseqFeatureCount <- function(dir.process, consensusToCount, sampleids, nthreads){
    # GET REGION TO COUNT ---
    regionsToCount <- data.frame(GeneID = paste(
                                            GenomicRanges::seqnames(consensusToCount), 
                                            GenomicRanges::start(consensusToCount), 
                                            GenomicRanges::end(consensusToCount), sep = "_"), 
                                Chr = GenomicRanges::seqnames(consensusToCount), 
                                Start = GenomicRanges::start(consensusToCount), 
                                End = GenomicRanges::end(consensusToCount), 
                                Strand = GenomicRanges::strand(consensusToCount))

    # GET FILE FILES TO COUNT ---
    bamfile_suffix <- ".bowtie.sorted.nodup.shifted.bam"
    files.bam <- sapply(sampleids, function(x) {
                        paste(dir.process, "/", sprintf("%s/alignment/%s", x, x), bamfile_suffix, sep="")
                    })

    # COUNT READS
    fcResults <- Rsubread::featureCounts(
                                files=files.bam, 
                                annot.ext=regionsToCount, 
                                isPairedEnd=TRUE, 
                                countMultiMappingReads=FALSE, 
                                maxFragLength=100,
                                autosort=TRUE,
                                nthreads=nthreads,
                                verbose=FALSE)
  
    # COLLECT THE RAW COUNTS ---
    counts <- fcResults$counts
    colnames(counts) <- sampleids
    counts <- cbind.data.frame(feature_id=rownames(counts), counts)
    rownames(counts) <- NULL

    return(counts)
}


### FUNCTION: get.NormFeatureExpression() ---
get.NormFeatureExpression <- function(dds){
    # LOG TRANSFORM DESEQ dds OBJECT ---
    if(nrow(dds@colData) < 40){
        rds <- DESeq2::rlog(dds)
    }else{
        rds <- DESeq2::vst(dds)
    }
    
    mat.rds <- SummarizedExperiment::assay(rds)

    return(mat.rds)
}



#' featurefilter: A function for filtering features
#'
#' This function is to filter features based on variance. Depending on the data different
#' metrics will be more appropiate, simple variance is included if variance does not tend to
#' increase with the mean. There is also the median absolute deviation which is a more robust
#' metric than variance, this is preferable. The coefficient of variation (A) or its second
#' order derivative (A2) (Kvalseth, 2017) are also included which standardise the standard
#' deviation with respect to the mean. It is best to manually examine the mean-variance relationship of 
#' the data, for example, using the results from this function together with the qplot function 
#' from ggplot2.
#'
#' @param mydata Data frame: should have samples as columns and rows as features
#' @param percentile Numerical value: the top X percent most variable features should be kept
#' @param method Character vector: variance (var), coefficient of variation (A), second order A (A2), median absolute deviation (MAD)
#' @param topN Numerical value: the number of most variable features to display
#'
#' @return A list, containing: 
#' 1) filtered data
#' 2) statistics for each feature order according to the defined filtering metric
#' 
#' @references 
#' Kvålseth, Tarald O. "Coefficient of variation: the second-order alternative." Journal of Applied Statistics 44.3 (2017): 402-415.
#' 
#' @export
#'
#' @examples
#' filtered <- featurefilter(mydata,percentile=10)
#' NOTE: THIS FUNCTION IS EXTRACTED FROM M3C PACKAGE https://bioconductor.org/packages/release/bioc/html/M3C.html
### FUNCTION: featurefilter() ---
featurefilter <- function(mydata,percentile=10,method='MAD',topN=20){
  
    message('***feature filter function***')
  
    message(paste('extracting the most variable: '),percentile,' percent')
    message(paste('features to start with:',nrow(mydata)))
  
    # percentile (convert to decimal below)
    percentile <- 1-(percentile/100)
  
    if (method == 'A'){
        message('performing calculations for co efficient of variation/A')
        
        # calculate mean and variance
        u <- rowMeans(mydata)
        sigma <- apply(mydata,1,sd)
        vars <- sigma^2

        # calc coefficient of variation for all rows (features)
        CV <- sigma/u
        CV[is.na(CV)] <- 0
        A <- CV

        # get features with CV in the given percentile
        CVthresh <- quantile(CV, percentile, na.rm = TRUE) 

    }else if (method == 'A2'){
        message('performing calculations for second order co efficient of variation/A2')
    
        # calculations
        u <- rowMeans(mydata)
        sigma <- apply(mydata,1,sd)
        vars <- sigma^2
        A <- sigma/u
        A[is.na(A)] <- 0
        AA <- A^2

        # get second order co efficient of variation
        A2 <- sqrt((AA/(AA+1)))
        CV <- A2

        # get features with CV in the given percentile
        CVthresh <- quantile(CV, percentile, na.rm = TRUE)

    }else if (method == 'var'){
        message('performing calculations for variance')

        u <- rowMeans(mydata)
        sigma <- apply(mydata,1,sd)
        vars <- sigma^2
        CV <- vars
        CVthresh <- quantile(CV, percentile, na.rm = TRUE)

    }else if (method == 'MAD'){
        message('performing calculations for median absolute deviation')

        u <- rowMeans(mydata)
        MAD <- apply(mydata,1,mad)
        sigma <- apply(mydata,1,sd)
        vars <- sigma^2
        CV <- MAD
        CVthresh <- quantile(CV, percentile, na.rm = TRUE)
    }
  
    ## filter data
    names <- names(CV)[CV>=as.numeric(CVthresh)]
    filtered_data <- subset(mydata, row.names(mydata) %in% names)
  
    # make data frame of results
    if (method == 'A'){
        test <- data.frame('feature'=row.names(mydata),'mean'=u,'var'=vars,'sd'=sigma,'A'=A)
        test <- test[order(-test[,5]), ]
        message('printing topN most variable features with statistics...')
        print(head(test,topN))
    }else if (method == 'A2'){
        test <- data.frame('feature'=row.names(mydata),'mean'=u,'var'=vars,'sd'=sigma,'A'=A,'A2'=A2)
        test <- test[order(-test[,6]), ]
        message('printing topN most variable features with statistics...')
        print(head(test,topN))
    }else if (method == 'var'){
        test <- data.frame('feature'=row.names(mydata),'mean'=u,'var'=vars,'sd'=sigma)
        test <- test[order(-test[,3]), ]
        message('printing topN most variable features with statistics...')
        print(head(test,topN))
    }else if (method == 'MAD'){
        test <- data.frame('feature'=row.names(mydata),'mean'=u,'var'=vars,'sd'=sigma,'MAD'=MAD)
        test <- test[order(-test[,5]), ]
        message('printing topN most variable features with statistics...')
        print(head(test,topN))
    }
  
    ## make results list
    mylist <- list('filtered_data'=filtered_data,'statistics'=test)
  
    ## message
    message(paste('features remaining:',nrow(filtered_data)))
  
    ##
    return(mylist)
}



### FUNCTION: gmtPathways() ---
#' Returns a list of pathways from a GMT file.
#' @param gmt.file Path to a GMT file.
#' @return A list of vectors with gene sets.
gmtPathways <- function(file.gmt){
    pathwayLines <- strsplit(readLines(file.gmt), "\t")
    pathways <- lapply(pathwayLines, tail, -2)
    names(pathways) <- sapply(pathwayLines, head, 1)
    return(pathways)
}


### FUNCTION: RankNorm() ---
# SOURCE: 'RNOmni' R-package https://cran.r-project.org/web/packages/RNOmni
# Purpose: Rank normal transform
# Updated: 19/10/11

#' Rank-Normalize
#'
#' Applies the rank-based inverse normal transform (INT) to a numeric vector.
#' The INT can be broken down into a two-step procedure. In the first, the
#' observations are transformed onto the probability scale using the empirical
#' cumulative distribution function (ECDF). In the second, the observations are
#' transformed onto the real line, as Z-scores, using the probit function.
#' @param u Numeric vector.
#' @param k Offset. Defaults to (3/8), correspond to the Blom transform.
#' @return Numeric vector of rank normalized measurements.
#' 
#' @importFrom stats qnorm
#' @export
#' @seealso 
#' \itemize{
#'   \item Direct INT test \code{\link{DINT}}.
#'   \item Indirect INT test \code{\link{IINT}}.
#'   \item Omnibus INT test \code{\link{OINT}}.
#' }
#'
#' @examples
#' # Draw from chi-1 distribution
#' y <- rchisq(n = 1e3, df = 1)
#' # Rank normalize
#' z <- RankNorm(y)
#' # Plot density of transformed measurement
#' plot(density(z))
#'
RankNorm <- function(u, k = 0.375){
    # Input checks. 
    if(!is.vector(u)){
        stop("A numeric vector is expected for u.")
    }

    if((k < 0) || (k > 0.5)){
        stop("Select the offset within the interval (0,0.5).")
    }
    
    if(sum(is.na(u)) > 0){
        stop("Please exclude observations with missing measurements.")
    }

    # Observations.
    n <- length(u)
  
    # Ranks.
    r <- rank(u)
  
    # Apply transformation.
    out <- qnorm((r - k) / (n - 2 * k + 1))

    return(out)
}



#########################################
#   --------  CONFUSION MATRIX -------- #  
#                                       #
#                   OBSERVED CLASS      #
#                       1       0       #
#                   -----------------   #
#   PREDICTED   1   |   TP  |   FP  |   #
#   CLASS       0   |   FN  |   TN  |   #
#                   -----------------   #
######################################### 

### FUNCTION: computePerformance() ---
computePerformance <- function(mat){
    # GET CLASS METRIC ---
    TP <- mat[1,1]
	FP <- mat[1,2]
    FN <- mat[2,1]
	TN <- mat[2,2]

    # COMPUTE PERFORMANCE METRICS ---
    sn <- TP/(TP+FN) # RECALL (TPR)
    sp <- TN/(FP+TN) # (TNR)
    ppv <- TP/(TP+FP) # PRECISION
    npv <- TN/(FN+TN)

    accuracy <- (TP+TN)/(TP+FP+FN+TN)
    f1_score <- 2*TP / ( (2*TP) + FP + FN )

	lrp <- sn/(1-sp) # likelihood ratio positive
	lrn <- (1-sn)/sp # likelihood ratio negative

	auc <- 0.5 * (sn + ppv)
					
	fpr <- FP / (FP+TN) # (1-specificity)
	fnr <- FN / (FN+TP)	

	# Matthews Correlation Coefficient (MCC) ------
    #mcc <- suppressWarnings( ( (TP*TN) - (FP*FN) ) / sqrt( (TP+FP)*(TP+FN)*(FP+TN)*(TN+FN) ) )
    mcc <- suppressWarnings( ( (TP*TN) - (FP*FN) ) / ( sqrt(TP+FP) * sqrt(TP+FN) * sqrt(FP+TN) * sqrt(TN+FN) ) )
	if(is.na(mcc)){
		mcc <- 0
	}

    # PREPARE OUTPUT DATA ---
    list.metric <- list()

    list.metric[[ "metric" ]][[ "TP" ]] <- TP
    list.metric[[ "metric" ]][[ "FP" ]] <- FP
    list.metric[[ "metric" ]][[ "FN" ]] <- FN
    list.metric[[ "metric" ]][[ "TN" ]] <- TN

    list.metric[[ "sensitivity" ]] <- sn
    list.metric[[ "specificity" ]] <- sp
    list.metric[[ "ppv" ]] <- ppv
    list.metric[[ "npv" ]] <- npv
    list.metric[[ "accuracy" ]] <- accuracy
    list.metric[[ "f1_score" ]] <- f1_score
    list.metric[[ "lrp" ]] <- lrp
    list.metric[[ "lrn" ]] <- lrn
    list.metric[[ "auc" ]] <- auc
    list.metric[[ "fpr" ]] <- fpr
    list.metric[[ "fnr" ]] <- fnr
    list.metric[[ "mcc" ]] <- mcc

    return(list.metric)
}



### FUNCTION: get.normalizeQuantile() ---
get.normalizeQuantile <- function(dat){
	dQNorm <- data.frame(aroma.light::normalizeQuantile(as.matrix(dat)))
    colnames(dQNorm) <- colnames(dat)
	return(dQNorm)
}	


#### FUNCTION: get_over_representation_analysis() ---
get_over_representation_analysis <- function(list.genesets, genes.queryset, genes.refset, p.threshold){
    # CREATE RESULT DATAFRAME ---
	dat.result <- data.frame(Category=as.character(names(list.genesets)))

    # COMPUTE ONE-TAIL FISHER EXACT TEST ---
	for(i in 1:length(list.genesets)){
		genes.genesets <- sort(unique(unlist( list.genesets[[i]] )), decreasing=FALSE)
		genes.interest <- intersect(genes.queryset, genes.genesets)
  
		yy <- length(intersect(genes.genesets, genes.interest))
		yn <- length(intersect(genes.genesets, setdiff(genes.refset, genes.interest)))
		ny <- length(intersect(setdiff(genes.refset, genes.genesets), genes.interest))
		nn <- length(intersect(setdiff(genes.refset,genes.interest), setdiff(genes.refset, genes.genesets)))
  
		fisherRes <- fisher.test(rbind(c(yy,yn),c(ny,nn)), alternative="greater")
		dat.result$pvalue[i] <- fisherRes$p.value
		dat.result$fdr[i] <- NA

		dat.result$overlap.percent[i] <- round(length(genes.interest)/length(genes.genesets) * 100, digit=2)
		dat.result$overlap.genes[i] <- paste(genes.interest, collapse=":")
	}

	# MULTIPLE TEST CORRECTION ---
	dat.result$fdr <- p.adjust(dat.result$pvalue, method="BH")
	
    # TRIM RESULTS ---
	dat.result <- dat.result[order(dat.result$pvalue, decreasing=F),]
	dat.result <- subset(dat.result, dat.result$fdr <= p.threshold)
  
    # FORMAT PVALUES ---
	#dat.result$pvalue <- format(dat.result$pvalue, scientific = TRUE, digits = 4)
	#dat.result$fdr <- format(dat.result$fdr, scientific = TRUE, digits = 4)
	
	return(dat.result)
}  


### FUNCTION: get.wilcox.rank.test() WILCOXN RANK SUM TEST ----
get.wilcox.rank.test <- function(dat, class1, class2){
	cat("COMPUTE WILCOX RANK TEST ...", "\n", sep="")
	#dat <- dat[-which(is.na(dat)),]
	
	#Get class indices
	index1 <- which(colnames(dat) %in% class1)
	index2 <- which(colnames(dat) %in% class2)
	
	# Check if the data points are constant ----
	exp <- dat[,c(index1, index2)]
	len <- apply(exp, 1, function(x) length(unique(as.numeric(x))))
	del.index <- which(len == 1)

	if(length(del.index) != 0){
		dat <- dat[-del.index,]
	}
	
	# Compute Median and StdDev for the Sample Class
	median1 <- apply(as.matrix(dat)[,index1], 1, median, na.rm=TRUE)
	median2 <- apply(as.matrix(dat)[,index2], 1, median, na.rm=TRUE)
	cat("MEDIAN COMPUTED ...", "\n", sep="")
	
	stdev1 <- apply(as.matrix(dat)[,index1], 1, sd, na.rm=TRUE)
	stdev2 <- apply(as.matrix(dat)[,index2], 1, sd, na.rm=TRUE)
	cat("STANDARD DEVIATION COMPUTED ...", "\n", sep="")
	
	#Compute FoldChange
	foldchange <- ifelse(median2 > median1, 2^median2/2^median1, -2^median1/2^median2)
	cat("FOLD CHANGE COMPUTED ...", "\n", sep="")
	
	cat("PERFORMING WILCOX RANK TEST ...", "\n", sep="")
	#Wilcox.rank.test
	func.wilcoxtest <- function(x, index1, index2){
		wilcoxtest <- wilcox.test(x[index1], x[index2], alternative="two.sided", na.action=na.omit)
		return(wilcoxtest)
	}

	wilcoxtest.list <- apply(as.matrix(dat), 1, function(x) func.wilcoxtest(x, index1, index2))
	cat("WILCOX RANK TEST COMPUTED ...", "\n", sep="")
	
	pvalue <- do.call(c,lapply(wilcoxtest.list,function(x) x$p.value))

	index0 <- which(pvalue == 0)
	if(length(index0) != 0){
		pvalue[index0] <- .Machine$double.eps #smallest value
	}

	fdr <- p.adjust(pvalue, method = "BH", n = length(pvalue))
	cat("P-VALUE COMPUTED ...", "\n", sep="")

	dat.summary <- data.frame(Gene=rownames(dat), 
								MedianA=median1, StdevA=stdev1, 
								MedianB=median2, StdevB=stdev2, 
								FoldChange=foldchange, 
								pvalue=pvalue, fdr=fdr)
	dat.summary$Gene <- as.character(dat.summary$Gene)
	dat.summary <- dat.summary[order(dat.summary$pvalue, decreasing=FALSE),]
    rownames(dat.summary) <- NULL
	cat("DATA SUMMARY COMPUTED ...", "\n", sep="")
		
	return(dat.summary)
}

##################################################################################
cat(format(Sys.time(), "%a %b %d %X %Y"), "DONE!", "\n", sep=" ")
##################################################################################
