###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### LOAD LIBRARIES ---
cat(format(Sys.time(), "%a %b %d %X %Y"), "LOADING PACKAGE: Sushi", "\n", sep=" ")
suppressPackageStartupMessages(require("Sushi", quietly=TRUE, warn.conflicts=FALSE))


##################################################################################
cat(format(Sys.time(), "%a %b %d %X %Y"), "LOADING FUNCTIONS ...", "\n", sep=" ")
##################################################################################

### FUNCTION: get.overlap() ---
# THIS FUNCTION RETURNS WHETHER TWO PAIRS OF LOCI OVERLAP
# y = START COORD
# z = END COORD
get.overlap <- function(x1, x2, y, z){
    stopifnot(length(x1)==length(x2))

    if(length(x1)==1 && !is.na(x1) && !is.na(x2) && x1>x2) {
        x3 <- x2
        x2 <- x1
        x1 <- x3
    }else if(length(x1)>1) {
        rowswrong <- x1>x2
        rowswrong[is.na(rowswrong)] <- F
        x3 <- x2[rowswrong]
        x2[rowswrong] <- x1[rowswrong]
        x1[rowswrong] <- x3
    }
  
    stopifnot(y<=z)

    x1in <- x1>=y & x1<=z
    x2in <- x2>=y & x2<=z
    xencompasses <- x1<y & x2>z
    
    return(x1in | x2in | xencompasses)
}


### FUNCTION: plot_gene_track() ---
plot_gene_track <- function(hg38_models, chrom_name, chrom_start, chrom_end){
  # GET EXON BED BY QUERY GENE ---
  exonrows <- hg38_models$Chrom==chrom_name & get.overlap(x1=hg38_models$Start, x2=hg38_models$End, y=chrom_start, z=chrom_end)
  hg_tmp <- hg38_models[exonrows,]

  if( dim(hg_tmp)[1] > 0 ){
    pg <- Sushi::plotGenes(geneinfo=hg_tmp, chrom=chrom_name, chromstart=chrom_start, chromend=chrom_end,
                    plotgenetype="box", type=rep("exon", dim(hg_tmp)[1]),
                    bheight=0.01, 
                    labeltext=TRUE, 
                    packrow=TRUE,
                    col="#000000",
                    fontsize=0.5, labeloffset = 0.5,
                    fonttype=3)

  }else{
    plot(-1,-1, xlim=c(0,10), ylim=c(0,10), xlab="", ylab="", axes=FALSE)   
  }
}

### FUNCTION: plot_mono_gene_track() ---
plot_mono_gene_track <- function(hg38_models, gene, chrom_name, chrom_start, chrom_end){
  # EXTRACT EXON MODEL BY GENE ---
  hg_tmp <- subset(hg38_models, hg38_models$Gene == gene)

  # PLOT ---
  if( nrow(hg_tmp) > 0 ){
    pg <- Sushi::plotGenes(geneinfo=hg_tmp, chrom=chrom_name, chromstart=chrom_start, chromend=chrom_end,
                    plotgenetype="box", type=rep("exon", nrow(hg_tmp)),
                    bheight=0.15, 
                    labeltext=TRUE, 
                    packrow=FALSE,
                    fontsize=1, labeloffset = 0.5,
                    fonttype=3)

  }else{
    plot(-1,-1, xlim=c(0,10), ylim=c(0,10), xlab="", ylab="", axes=FALSE)   
  }
}



### FUNCTION: plot_coverage_track() ---
plot_coverage_track <- function(gr_bedgraph, chrom_name, chrom_start, chrom_end, color="", label_text="", range=NULL, textyPos=NULL){
  if(color == ""){
    color <- Sushi::SushiColors(2)(2)[1]
  }

  # BEDGRAPH DATAFRAME ---    
  data <- as.data.frame(gr_bedgraph)
  data <- subset(data, select=c("seqnames","start","end","score"))
  colnames(data) <- c("chrom","start","end","value")

  # CREATE UPPER BOUND FOR SCORE ---
  #if(max(data$score) > 100){
  #  data$score[which(data$score >= 100)] <- 100
  #}

  # PLOT --
  Sushi::plotBedgraph(signal=data, chrom=chrom_name, chromstart=chrom_start, chromend=chrom_end, lwd=1, range=range, transparency=.50, color=color, flip=FALSE)

  axis(side=2, las=2, tcl=.1, labels=TRUE, cex.axis=0.5, lwd=0.5)
  #for( x in seq(from=0,to=2500,by=500)){
  #for( x in seq(from=1,to=10,by=1)){
    #abline( x,0, col="#0000ff33")   
  #}
  abline( 0,0, col="black")   
  
  if( length(label_text)>0 ){
    maxy <- 0.9
    if( dim(data)[1]>0 ){
      dvals <- data$value[data$start>=chrom_start & data$end<=chrom_end & data$chrom==chrom_name]
      if( length(dvals)>0 ){
        maxy <- max(dvals, na.rm=TRUE) * 0.9
      }
    }
    
    if(is.null(textyPos)){
      textyPos = maxy
    }
    
    text( chrom_start + (chrom_end-chrom_start)/50, textyPos, 
          adj=0, label_text, font=1, col="#000000", cex=0.5)
  }  
}




### FUNCTION: plotBed_track() ---
plotBed_track <- function(gr_bed_peak, chrom_name, chrom_start, chrom_end, label_text="", range=NULL, textyPos=NULL){
  # PEAK BED DATAFRAME ---    
  if(length(gr_bed_peak) != 0){
    data <- as.data.frame(gr_bed_peak)
    data <- subset(data, select=c("seqnames","start","end"))
    colnames(data) <- c("chrom","start","end")

    # PLOT --
    Sushi::plotBed(beddata=data, chrom=chrom_name, chromstart=chrom_start, chromend=chrom_end,
            type = "region", height=1, wiggle = 0, color=Sushi::SushiColors(2)(2)[2])

  }else if(length(gr_bed_peak) == 0){
    gr_plotrange <- GenomicRanges::GRanges(seqnames = Rle(chrom_name), ranges = IRanges(start=chrom_start, end=chrom_end))
    gr_plotrange_bins <- GenomicRanges::tile(x=gr_plotrange, width=10)
    data <- as.data.frame(gr_plotrange_bins[[1]][2])[,1:3]
    colnames(data) <- c("chrom","start","end")

    # PLOT --
    Sushi::plotBed(beddata=data, chrom=chrom_name, chromstart=chrom_start, chromend=chrom_end,
            type = "region", height=1, wiggle = 0.02, color=Sushi::opaque(color = Sushi::SushiColors(2)(2)[1], transparency = 0))
  }



  # PLOT --
  #Sushi::plotBed(beddata=data, chrom=chrom_name, chromstart=chrom_start, chromend=chrom_end,
  #          type = "region", height=2, wiggle = 0.02, color=Sushi::SushiColors(2)(2)[1])


  #axis(side=2,las=2,labels=FALSE,tick=FALSE )
  #for( x in seq(from=1,to=10,by=1)){
  #  abline( x,0, col="#0000ff33")   
  #}
  abline( 0,0, col="black")   
  


  if( length(label_text)>0 ){
    maxy <- 0.9
    if( dim(data)[1]>0 ){
      dvals = data$value[data$start>=chrom_start & data$end<=chrom_end & data$chrom==chrom_name]
      if( length(dvals)>0 ){
        maxy <- max(dvals, na.rm=TRUE) * 0.9
      }
    }
    
    if(is.null(textyPos)){
      textyPos = maxy
    }
    
    text( chrom_start + (chrom_end-chrom_start)/50, textyPos, 
          adj=0, label_text, font=2, col="darkblue")
  }  
}



#########################################################################################################


### FUNCTION: loadBedGraph() ---
loadBedGraph <- function(files.bedgraph, n_cores){
    list.bedgraph <- list()
    for(k in 1:length(files.bedgraph)){
        subtype <- names(files.bedgraph)[k]
        file.bedgraph <- files.bedgraph[subtype]

        # LOAD BEDGRAPH ---
        cat(format(Sys.time(), "%a %b %d %X %Y"), "LOADING BEDGRAPH:", subtype, "\n", sep=" ")
            dt_bedgraph <- data.table::fread(file.bedgraph, data.table=TRUE, nThread=n_cores)
            colnames(dt_bedgraph) <- c("chr","start","end","score")

        # CONVERT TO GRANGES OBJECT ---
        cat(format(Sys.time(), "%a %b %d %X %Y"), "CONVERTING BEDGRAPH TO GRANGE OBJECT", "\n", sep=" ")
            list.bedgraph[[subtype]] <- GenomicRanges::makeGRangesFromDataFrame(dt_bedgraph, keep.extra.columns = TRUE)
        cat(format(Sys.time(), "%a %b %d %X %Y"), "DONE!", "\n", sep=" ")
    }

    return(list.bedgraph)
}

### FUNCTION: loadPeaks() ---
loadPeaks <- function(files.peaks, n_cores){
    list.peaks <- list()
    for(k in 1:length(files.peaks)){
        subtype <- names(files.peaks)[k]
        file.bed <- files.peaks[subtype]

        # LOAD BEDGRAPH ---
        cat(format(Sys.time(), "%a %b %d %X %Y"), "LOADING PEAKS:", subtype, "\n", sep=" ")
            dt_bed <- data.table::fread(file.bed, data.table=TRUE, nThread=n_cores)
            colnames(dt_bed) <- c("chr","start","end")

        # CONVERT TO GRANGES OBJECT ---
        cat(format(Sys.time(), "%a %b %d %X %Y"), "CONVERTING PEAKS TO GRANGE OBJECT", "\n", sep=" ")
            list.peaks[[subtype]] <- GenomicRanges::makeGRangesFromDataFrame(dt_bed, keep.extra.columns = FALSE)
        cat(format(Sys.time(), "%a %b %d %X %Y"), "DONE!", "\n", sep=" ")
    }

    return(list.peaks)
}

### FUNCTION: loadDiffPeaks() ---
loadDiffPeaks <- function(file.bed){
    dat_diffpeaks <- data.table::fread(file=file.bed, sep="\t", header=FALSE, nThread=60, data.table=FALSE, verbose=FALSE)
    colnames(dat_diffpeaks) <- c("chr","start","end","gene")
    gr_diffpeaks <- GenomicRanges::makeGRangesFromDataFrame(df=dat_diffpeaks, keep.extra.columns=TRUE)
    return(gr_diffpeaks)
}


### FUNCTION: subsetPeaksByRegion() ---
subsetPeaksByRegion <- function(list.peaks, chrom_name, chrom_start, chrom_end){
    gr_region <- GenomicRanges::GRanges(seqnames = Rle(chrom_name), ranges = IRanges(start=chrom_start, end=chrom_end))

    sampleids <- names(list.peaks)
    list.peaks.region <- list()
    for(k in 1:length(list.peaks)){
        sampleid <- sampleids[k]
        gr_peak <- GenomicRanges::makeGRangesFromDataFrame(list.peaks[[sampleid]], keep.extra.columns = TRUE)
        
        # SUBSET PEAKS BY PLOT REGION ---
        gr_peak_region <- IRanges::subsetByOverlaps(gr_peak, gr_region)
        list.peaks.region[[sampleid]] <- gr_peak_region

        cat(format(Sys.time(), "%a %b %d %X %Y"), "PROCESSED:", sampleid, "\n", sep=" ")
    }

    return(list.peaks.region)
}


### FUNCTION: getCoverageByPeakRegion() ---
getCoverageByPeakRegion <- function(list.bedgraph, list.peaks){
    sampleids <- names(list.peaks)
    list.bedgraph.region <- list()
    for(k in 1:length(list.peaks)){
        sampleid <- sampleids[k]

        gr_bedgraph <- list.bedgraph[[sampleid]]
        gr_peaks <- list.peaks[[sampleid]]

        # SUBSET BEDGRAPH BY PLOT REGION ---
        gr_coverage_region <- IRanges::subsetByOverlaps(gr_bedgraph, gr_peaks)
        list.bedgraph.region[[sampleid]] <- gr_coverage_region

        cat(format(Sys.time(), "%a %b %d %X %Y"), "PROCESSED:", sampleid, "\n", sep=" ")
    }
    return(list.bedgraph.region)
}


### FUNCTION: plotBed_track_modified() ---
plotBed_track_modified <- function(gr_bed_peak, chrom_name, chrom_start, chrom_end, color=Sushi::SushiColors(2)(2)[2], label_text="", range=NULL, textyPos=NULL){
  # PEAK BED DATAFRAME ---    
  if(length(gr_bed_peak) != 0){
    data <- as.data.frame(gr_bed_peak)
    data <- subset(data, select=c("seqnames","start","end"))
    colnames(data) <- c("chrom","start","end")

    # PLOT --
    Sushi::plotBed(beddata=data, chrom=chrom_name, chromstart=chrom_start, chromend=chrom_end,
            type = "region", height=0.5, wiggle=0, color=color)

  }else if(length(gr_bed_peak) == 0){
    gr_plotrange <- GenomicRanges::GRanges(seqnames = Rle(chrom_name), ranges = IRanges(start=chrom_start, end=chrom_end))
    gr_plotrange_bins <- GenomicRanges::tile(x=gr_plotrange, width=10)
    data <- as.data.frame(gr_plotrange_bins[[1]][2])[,1:3]
    colnames(data) <- c("chrom","start","end")

    # PLOT --
    Sushi::plotBed(beddata=data, chrom=chrom_name, chromstart=chrom_start, chromend=chrom_end,
            type = "region", height=1, wiggle = 0, color=Sushi::opaque(color = Sushi::SushiColors(2)(2)[1], transparency = 0))
  }

  #axis(side=2,las=2,labels=FALSE,tick=FALSE )
  #for( x in seq(from=1,to=10,by=1)){
  #  abline( x,0, col="#0000ff33")   
  #}
  abline( 0,0, col="black")   
  
  if( length(label_text)>0 ){
    maxy <- 0.9
    if( dim(data)[1]>0 ){
      dvals = data$value[data$start>=chrom_start & data$end<=chrom_end & data$chrom==chrom_name]
      if( length(dvals)>0 ){
        maxy <- max(dvals, na.rm=TRUE) * 0.9
      }
    }
    
    if(is.null(textyPos)){
      textyPos = maxy
    }
    
    text( chrom_start + (chrom_end-chrom_start)/50, textyPos, 
          adj=0, label_text, font=0.5, col="#000000", cex=0.5)
  }  
}

### FUNCTION: plot_links() ---
plot_links <- function(bedped, plot_chr, plot_region_start, plot_region_end, label_text, textyPos=NULL){
    # PLOT: LINKS ----
    Sushi::plotBedpe(bedpedata=bedped, chrom=plot_chr, chromstart=plot_region_start, chromend=plot_region_end,
            heights=bedped$score, plottype="loops", color="#000000", lwd=0.5, flip=TRUE)

    abline( 0,0, col="black")  
    axis( side=2, las=2, tcl=0.1, labels=-1*pretty(par("yaxp")[c(1,2)]), at=pretty(par("yaxp")[c(1,2)]), lwd=0.5, lwd.ticks=0.5, tick=TRUE)
    text( plot_region_start + (plot_region_end-plot_region_start)/50, textyPos, adj=0, labels=label_text, font=0.2, cex=0.2, col="#000000")
}


### FUNCTION: getPlot() ---
getPlot <- function(annot.exon, list.bedgraph.region, gr_diffpeaks, gr_h3k27ac, bedped=NULL, plot_chr, plot_region_start, plot_region_end, subtypes, cp_subtypes, file.plot, height=4.5, width=2, range=c(0,1.5), textyPos=0.7){
    pdf(file.plot, height=height, width=width)
        # MANAGE PLOT AREA ---
        nrow <- 11
        heights <- NULL

        if(is.null(heights)){
            heights <- rep(1, nrow)
            heights[1] <- 2
            heights[ length(heights) ] <- 2
        }

        # PLOT START ----
        layout(matrix(1:nrow,nrow,1), heights=heights)
        par(mar=c(0,3,1,1))

        # PLOT GENE MODEL ---
        plot_gene_track(hg38_models=annot.exon, chrom_name=plot_chr, chrom_start=plot_region_start, chrom_end=plot_region_end)

        par(mar=c(1,3,1,1))
        # PLOT DIFFPEAKS TRACK ---
        plotBed_track_modified(gr_bed_peak=gr_diffpeaks, chrom_name=plot_chr, chrom_start=plot_region_start, chrom_end=plot_region_end, color="#000000", label_text="DiffPeaks", range=NULL, textyPos=NULL)

        par(mar=c(1,3,1,1))
        # PLOT ATAC-SEQ PEAKS ---
        gr_bedgraph_region <- list.bedgraph.region[[subtypes[1]]]
        plot_coverage_track(gr_bedgraph_region, chrom_name=plot_chr, chrom_start=plot_region_start, chrom_end=plot_region_end, color=cp_subtypes[1], label_text=subtypes[1], range=range, textyPos=textyPos)

        gr_bedgraph_region <- list.bedgraph.region[[subtypes[2]]]
        plot_coverage_track(gr_bedgraph_region, chrom_name=plot_chr, chrom_start=plot_region_start, chrom_end=plot_region_end, color=cp_subtypes[2], label_text=subtypes[2], range=range, textyPos=textyPos)

        gr_bedgraph_region <- list.bedgraph.region[[subtypes[3]]]
        plot_coverage_track(gr_bedgraph_region, chrom_name=plot_chr, chrom_start=plot_region_start, chrom_end=plot_region_end, color=cp_subtypes[3], label_text=subtypes[3], range=range, textyPos=textyPos)

        gr_bedgraph_region <- list.bedgraph.region[[subtypes[4]]]
        plot_coverage_track(gr_bedgraph_region, chrom_name=plot_chr, chrom_start=plot_region_start, chrom_end=plot_region_end, color=cp_subtypes[4], label_text=subtypes[4], range=range, textyPos=textyPos)

        gr_bedgraph_region <- list.bedgraph.region[[subtypes[5]]]
        plot_coverage_track(gr_bedgraph_region, chrom_name=plot_chr, chrom_start=plot_region_start, chrom_end=plot_region_end, color=cp_subtypes[5], label_text=subtypes[5], range=range, textyPos=textyPos)

        # PLOT H3K27ac TRACK ---
        plotBed_track_modified(gr_bed_peak=gr_h3k27ac, chrom_name=plot_chr, chrom_start=plot_region_start, chrom_end=plot_region_end, color="#000000", label_text="H3K27ac", range=NULL, textyPos=NULL)

        # PLOT H3K27ac HIChIP TRACK ---
        #plot_links(bedped, plot_chr, plot_region_start, plot_region_end, label_text="", textyPos=NULL)

        par(mar=c(1,3,1,1))
        # PLOT GENOME LABEL ---
        Sushi::labelgenome(chrom=plot_chr, chromstart=plot_region_start, chromend=plot_region_end, n=10, chromcex=0.5, scalecex=0.5, scale="Mb")
    dev.off()
}


### FUNCTION: getPlotbyGene() ---
getPlotbyGene <- function(annot.gene, annot.exon, list.peaks, list.bedgraph, diffpeaks_bed, gr_h3k27ac, bedped=NULL, file.plot, subtypes, cp_subtypes, upstream=50000, downstream=10000, height, width, range, textyPos, gene){
    # GET GENE POSITION ---
    plot_chr <- annot.gene$Chrom[which(annot.gene$Gene == gene)] 
    gene_start <- annot.gene$Start[which(annot.gene$Gene == gene)] 
    gene_end <- annot.gene$End[which(annot.gene$Gene == gene)]
    plot_region_start <- gene_start - upstream   #50000 bp
    plot_region_end <- gene_end + downstream     #10000 bp

    # RANGE
    gr_range_limit <- GenomicRanges::GRanges(
                                seqnames = Rle(plot_chr),
                                ranges = IRanges(start=plot_region_start, end=plot_region_end),
                                strand = Rle(strand("+")))

    # EXTRACT PEAKS BY PLOT REGION ---
    list.peaks.region <- subsetPeaksByRegion(list.peaks, chrom_name=plot_chr, chrom_start=plot_region_start, chrom_end=plot_region_end)

    # EXTRACT COVERAGE BY PLOT PEAKS REGION ---
    list.bedgraph.region <- getCoverageByPeakRegion(list.bedgraph, list.peaks=list.peaks.region)

    # GET H3K27ac CHIP-SEQ PEAKS OVERLAP WITH PLOT GRANGE ---
    gr_h3k27ac <- IRanges::subsetByOverlaps(x=gr_h3k27ac, ranges=gr_range_limit, type="any")
    gr_h3k27ac <- gr_h3k27ac[unique(findOverlaps(gr_h3k27ac, type = "any", select = "first"))]

    # LOAD ATAC-SEQ DIFFPEAKS  OVERLAP WITH PLOT GRANGE ---
    gr_diffpeaks <- GenomicRanges::makeGRangesFromDataFrame(df=diffpeaks_bed, keep.extra.columns=TRUE)
    gr_diffpeaks <- IRanges::subsetByOverlaps(x=gr_diffpeaks, ranges=gr_range_limit, type="any")

    # PLOT ---
    getPlot(annot.exon, list.bedgraph.region, gr_diffpeaks, gr_h3k27ac, bedped=NULL, plot_chr, plot_region_start, plot_region_end, subtypes, cp_subtypes, file.plot, height, width, range=range, textyPos=textyPos)
}

