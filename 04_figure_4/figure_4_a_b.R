###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### DEFINE LIBRARIES ---
library("stringr")
library("data.table")
library("dplyr")
library("ggplot2")
library("gridExtra")

### DEFINE PATH ---
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

dir.footprints <- file.path(dir.reproduce_data, "tf_footprints/mcrpc_subtypes")
#dir.footprints <- file.path(dir.wrk, "analysis/12_tf_footprint/subtypes_AR/tobias")

dir.data <- tempdir() # REPLACE PATH
dir.scripts <- tempdir() # REPLACE PATH
dir.logs <- tempdir() # REPLACE PATH

### PARAMS ---
subtype <- "ARpNEn"
tfs <- c("AR","HOXB13")
motifs <- c("ExtendedSite_AR","MA0901.1_HOXB13")
names(motifs) <- tfs

### DEFINE FILES ---
file.rds_footprints <- file.path(dir.footprints, "data_footprints/tobias_mcrpc_subtypes_tf_footprints_featureids_bound_unbound.rds")
files.bigwig_corrected <- sapply(subtype, function(x) { paste(dir.footprints, "/", sprintf("%s/%s", x, x), "_corrected.bw", sep="") })

files.matrix_bound <- sapply(tfs, function(x) { paste(dir.data, "/matrix_bound_", sprintf("%s", x), ".tsv.gz", sep="") })
files.matrix_unbound <- sapply(tfs, function(x) { paste(dir.data, "/matrix_unbound_", sprintf("%s", x), ".tsv.gz", sep="") })



### FUNCTION: getBED() ---
getBED <- function(feature_ids){
    # GET BED FILE ---
    df.bed <- data.frame(chr=unlist(lapply(stringr::str_split(feature_ids, "_"), function(x) x[1])),
                        start=as.numeric(unlist(lapply(stringr::str_split(feature_ids, "_"), function(x) x[2]))),
                        end=as.numeric(unlist(lapply(stringr::str_split(feature_ids, "_"), function(x) x[3]))))

    # REMOVE DUPLICATES ---
    df.bed <- df.bed[!duplicated(df.bed),]

    return(df.bed)
}

### FUNCTION: getMotifWidth() ---
getMotifWidth <- function(df.bed){
    gr <- GenomicRanges::makeGRangesFromDataFrame(df.bed, keep.extra.columns=FALSE)
    motif_width <- GenomicRanges::width(gr)[1]
    return(motif_width)
}


### FUNCTION: getComuteMatrix() ---
getComuteMatrix <- function(file.bigwig, file.bed, file.output_matrix, subtype, motif_width, n_cores){
    # 'computeMatrix' FUNCTION FROM 'deepTools'
    # https://deeptools.readthedocs.io/en/latest/index.html

    cmd <- paste("computeMatrix scale-regions",
                "--scoreFileName", file.bigwig,
                "--regionsFileName", file.bed,
                "--outFileName", file.output_matrix,
                "--regionBodyLength", motif_width,
                "--beforeRegionStartLength", 60,
                "--afterRegionStartLength", 60,
                "--binSize", 1,
                "--samplesLabel", subtype,
                "--missingDataAsZero",
                "--verbose",
                "-p", n_cores,  
                sep=" ")
    
    return(cmd)
}



### FUNCTION: prepareData() ---
loadData <- function(file.matrix, group){
    # PARSE MOTIF_WIDTH FROM FIRST LINE MAETADATA ---
    line_1 <- data.table::fread(file=file.matrix, sep="\t", header=FALSE, nThread=55,  data.table=FALSE, verbose=FALSE, nrow=1)$V1
    line_1_sptxt <- stringr::str_split(line_1, ",")[[1]]
    line_1_body <- stringr::str_replace_all(stringr::str_replace_all(line_1_sptxt[str_detect(line_1_sptxt, "body")], "\\[", ""), "\\]", "")
    motif_width <- as.numeric(unlist(lapply(stringr::str_split(line_1_body, ":"), function(x) x[2])))

    # LOAD DATA ---
    dat <- data.table::fread(file=file.matrix, sep="\t", header=FALSE, nThread=55,  data.table=FALSE, verbose=FALSE, skip=1)

    # GET FEATURE ---
    feature_id <- dat$V4
    feature_id <- stringr::str_replace_all(stringr::str_replace_all(feature_id, ":", "_"), "-", "_")

    # GET MATRIX ---
    mat <- dat[,-c(1:6)]
    rownames(mat) <- feature_id
    colnames(mat) <- c(1:ncol(mat))
    #colnames(mat) <- c( paste("n", c(60:1), sep="_"), c(1:motif_width), paste("p", c(1:60), sep="_") )


    # SCALE THE SIGNAL BY MEAN OF THE THE LEFT FLANKING REGION ---
    #y_motif_start <- apply(mat, 1, function(x) as.numeric(x[61]))
    y_motif_start <- apply(mat, 1, function(x) mean(as.numeric(x[1:60])) )
    mat1 <- mat - y_motif_start

    # RESHAPE DATA ---
    df <- reshape2::melt(as.matrix(mat1))
    colnames(df) <- c("FeatureID","Position","Value")

    # GET MEAN PROFILE ---
    dm <- df %>% 
            dplyr::group_by(Position) %>%
            dplyr::summarize(Mean=mean(Value), 
                                SD=sd(Value),
                                Upper=Mean+SD,
                                Lower=Mean-SD)
    # ADD GROUP ---
    dm$Group <- group

    return(dm)
}

### FUNCTION: prepareData() ---
prepareData <- function(file.matrix_bound, file.matrix_unbound){
    # LOAD DATA ---
    dm_bound <- loadData(file.matrix=file.matrix_bound, group="Bound")
    dm_unbound <- loadData(file.matrix=file.matrix_unbound, group="Unbound")

    # MERGE DATA ---
    dm <- rbind(dm_bound, dm_unbound)

    # FACTORIZE ---
    dm$Group <- factor(dm$Group, levels=c("Bound","Unbound"))

    return(dm)
}


### FUNCTION: getPlot() ---
getPlot <- function(dm, gene, subtype){
    # DEFINE TITLE ---
    title.name <- paste(gene, subtype, sep="  ")

    # DEFINE COLOR ---
    cpalette <- c("#E31A1C", "#000000")
    #fpalette <- c("#D9D9D9", "#D9D9D9")

    # GET MARKED LINE POSITION ---
    y1 <- sort(unique(dm$Position))
    motif_start_pos <- y1[60+1]
    motif_end_pos <- y1[length(y1)-60]


    # PLOT ---
    p <- ggplot(data = dm, aes(x=Position, y=Mean)) +
            geom_vline(xintercept=motif_start_pos, color="gray70", linetype=4, alpha=0.5, size=0.25) +		
            geom_vline(xintercept=motif_end_pos, color="gray70", linetype=4, alpha=0.5, size=0.25) +	
            #geom_ribbon(aes(ymin=Lower, ymax=Upper, fill=Group), alpha=0.2) +
            #geom_line(aes(y=Upper, color=Group), size=0.25, alpha=0.5, linetype=4) +
            #geom_line(aes(y=Lower, color=Group), size=0.25, alpha=0.5, linetype=4) +            
            geom_line(aes(y=Mean, color=Group), alpha=0.8, size=0.4) +
            #geom_errorbar(data=dm, aes(x=Sample, ymin=Upper, ymax=Lower), width=0.2, size=0.5, color="#969696") +
            #geom_point(color="#E31A1C", size=1) +
            scale_x_continuous(breaks=c(1, motif_start_pos, motif_end_pos, max(y1)), labels=c("-60", "s", "e", "60")) +
            scale_y_continuous(breaks=c(-0.15,-0.10,-0.05,0,0.05,0.1), labels=c(-0.15,-0.10,-0.05,0,0.05,0.1)) +
            coord_cartesian(ylim=c(-0.15,0.1)) +
            scale_color_manual(values=cpalette) +
            scale_fill_manual(values=cpalette) +
            theme(
                axis.text.x = element_text(size = 5, color="#000000"),
                axis.text.y = element_text(size = 5, color="#000000"),
                axis.title = element_text(size = 5, color="#000000"),
                plot.title = element_text(size = 5, color="#000000", hjust=0.5),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.ticks = element_line(size=0.2, color="#000000"),	
                strip.text = element_text(size=5, color="#000000"),
                strip.background = element_rect(fill="#FFFFFF", color="#FFFFFF"),
                panel.background = element_rect(fill="#FFFFFF", color="#000000"),
                legend.text = element_text(size = 5, color="#000000"),
                legend.title = element_blank(),
                legend.key.size = unit(0.3, "cm"),
                legend.position = "bottom") +
            ylab("") +
            xlab("") + 
            ggtitle(title.name)

    return(p)
}


############################################################################################################

### LOAD FOOTPRINTS ---
list.rds_footprints <- readRDS(file=file.rds_footprints)



### GET FOOTPRINTS BED ---
list.bed_bound <- list.bed_unbound <-list()
list.files_bed_bound <- list.files_bed_unbound <- list()
list.motif_width_bound <- list.motif_width_unbound <-list()
for(tf in tfs){
    cat("GENERATING BED:", tf, "\n", sep=" ")

    motif_prefix <- motifs[tf]
    list.bed_bound[[tf]] <- getBED(feature_ids=unique( unlist(list.rds_footprints$bound[[motif_prefix]]) ))
    list.bed_unbound[[tf]] <- getBED(feature_ids=unique( unlist(list.rds_footprints$unbound[[motif_prefix]]) ))

    # GET MOTIF WITDH ---
    list.motif_width_bound[[tf]] <- getMotifWidth(df.bed=list.bed_bound[[tf]])
    list.motif_width_unbound[[tf]] <- getMotifWidth(df.bed=list.bed_unbound[[tf]])

    # DEFINE BED FILES ---
    list.files_bed_bound[[tf]] <- file.path(dir.scripts, paste("footprints_bound_", tf, ".bed", sep=""))
    list.files_bed_unbound[[tf]] <- file.path(dir.scripts, paste("footprints_unbound_", tf, ".bed", sep=""))

    # WRITE BED FILE ---
    write.table( list.bed_bound[[tf]], list.files_bed_bound[[tf]], sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
    write.table( list.bed_unbound[[tf]], list.files_bed_unbound[[tf]], sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

    cat("DONE:", tf, "\n", sep=" ")
}


### COMPUTE MATRIX ---  
for(tf in tfs){
    cat("COMPUTING MATRIX:", tf, "\n", sep=" ")

    cmd_bound_1 <- getComuteMatrix(file.bigwig=files.bigwig_corrected[subtype], 
                            file.bed=list.files_bed_bound[[tf]], 
                            file.output_matrix=files.matrix_bound[tf], 
                            subtype=subtype, 
                            motif_width=list.motif_width_bound[[tf]], 
                            n_cores=50)

    cmd_unbound_2 <- getComuteMatrix(file.bigwig=files.bigwig_corrected[subtype], 
                            file.bed=list.files_bed_unbound[[tf]], 
                            file.output_matrix=files.matrix_unbound[tf], 
                            subtype=subtype, 
                            motif_width=list.motif_width_bound[[tf]], 
                            n_cores=50)

    # EXECUTE COMMAND ---
    system(cmd_bound_1)
    system(cmd_unbound_2)

    cat("PROCESSED:", tf, "\n", sep=" ")
}




############################################################################################################

### GENERATE PLOT ---  
list.plot <- list()
for(tf in tfs){
    cat("START:", tf, "\n", sep=" ")

    dm <- prepareData(file.matrix_bound=files.matrix_bound[tf], 
                        file.matrix_unbound=files.matrix_unbound[tf])

    list.plot[[tf]] <- getPlot(dm, gene=tf, subtype)

    cat("PROCESSED:", tf, "\n", sep=" ")
}


### PLOT ---
file.plot <- file.path(dir.reproduce_fig, "figure_4_a_b.pdf")
pdf(file.plot, width=2, height=4)
        grid.arrange(list.plot[[ tfs[1] ]], list.plot[[ tfs[2] ]], nrow=2, ncol=1)
dev.off()
