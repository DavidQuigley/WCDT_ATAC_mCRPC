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
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/analysis/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

dir.bigwig <- file.path(dir.wrk, "secondary_data/wcdt_atacseq/bigwig") # REPLACE PATH
dir.matrix <- tempdir() # REPLACE PATH
dir.scripts <- tempdir() # REPLACE PATH
dir.logs <- tempdir() # REPLACE PATH

### DEFINE FILES ----
file.rds_atacseq_masterdata <- file.path(dir.reproduce_data, "wcdt_atacseq_masterdata.rds")
file.bed <- file.path(dir.wrk, "reference/annotation_data/gene_annotation/gencode_v28/genes_expressed_gencode_v28.bed")

### FUNCTION: get.SGEcmd()---
get.SGEcmd <- function(jobid, dir.logs){
    jobid <- jobid
    line01 <- "#!/bin/bash"
    line02 <- "#$ -pe smp 5"
    line03 <- "#$ -V"
    line04 <- "#$ -R y"
    line05 <- "#$ -l mem_free=10G"
    line06 <- paste("#$ -N", jobid, sep=" ")
    line07 <- paste("#$ -o", file.path(dir.logs, paste(jobid, ".out", sep="")), sep=" ")
    line08 <- paste("#$ -e", file.path(dir.logs, paste(jobid, ".err", sep="")), sep=" ")

    # COMPILE CMD ---
    cmd <- list(line01, line02, line03, line04, line05, line06, line07, line08)

    return(cmd)    
}

### FUNCTION: getComuteMatrix() ---
getComuteMatrix <- function(files.bigwig, file.bed, file.matrix, sampleids, n_cores=10){
    # 'computeMatrix' FUNCTION FROM 'deepTools'
    # https://deeptools.readthedocs.io/en/latest/index.html

    cmd <- paste("computeMatrix reference-point",
                "--scoreFileName", files.bigwig,
                "--regionsFileName", file.bed,
                "--outFileName", file.matrix,
                "--referencePoint", "TSS",
                "--beforeRegionStartLength", 3000,
                "--afterRegionStartLength", 3000,
                "--binSize", 100,
                "--samplesLabel", sampleids,
                "--missingDataAsZero",                
                "--verbose",
                "-p", n_cores,  
                sep=" ")
    
    return(cmd)
}


### LOAD ATAC-SEQ MASTER DATA ---
list.rds_atacseq_masterdata <- readRDS(file=file.rds_atacseq_masterdata)
sampleids <- list.rds_atacseq_masterdata$metadata_mcrpc$Sample_Name

### GET FILES ---
files.bigwig <- sapply(sampleids, function(x) { paste(dir.bigwig, "/", sprintf("%s", x), ".bw", sep="") })


### LOOP FOR EACH SAMPLE ---
list.sh <- list()
for(sampleid in sampleids){
    # GET FILES ---
    file.bw <- files.bigwig[sampleid]
    file.matrix <- file.path(dir.matrix, paste("computematrix_", sampleid, ".tsv.gz", sep="" ))

    # GET COMMANDS ---
    list.cmd_sge <- get.SGEcmd(jobid=paste("m", sampleid, sep="_"), dir.logs)
    cmd_cm <- getComuteMatrix(files.bigwig=file.bw, file.bed=file.bed, file.matrix=file.matrix, sampleids=sampleid, n_cores=55)

    # MERGE COMMANDS ---
    cmd <- c(unlist(list.cmd_sge), cmd_cm)

    # WRITE QSUB SCRIPT --
    file.qsub <- file.path(dir.scripts, paste("qsub_", sampleid, ".sh", sep=""))
    write.table(cmd, file.qsub, row.names=FALSE, col.names=FALSE, quote=FALSE)

    list.sh[[sampleid]] <- paste("qsub -q main", file.qsub, sep=" ")

    cat("PROCESSED:", sampleid, "\n", sep=" ")
}

### WRITE MAIN RUN SCRIPT ---
file.sh <- file.path(dir.scripts, "main.sh")
write.table(unlist(list.sh), file.sh, row.names=FALSE, col.names=FALSE, quote=FALSE)

cat("FILE GENERATED:", file.sh, "\n", sep=" ")






########################################################################################################################################################
# PLOT
########################################################################################################################################################
### FUNCTION: loadData() ---
loadData <- function(file.matrix, sampleid){
    # LOAD DATA ---
    dat <- data.table::fread(file=file.matrix, sep="\t", header=FALSE, nThread=55,  data.table=FALSE, verbose=FALSE, skip=1)

    # GET FEATURE ---
    feature_id <- dat$V4

    # GET MATRIX ---
    mat <- dat[,-c(1:6)]
    rownames(mat) <- feature_id
    colnames(mat) <- c(1:ncol(mat))

    # RESHAPE DATA ---
    df <- reshape2::melt(as.matrix(mat))
    colnames(df) <- c("FeatureID","Position","Value")

    # GET MEAN PROFILE ---
    dm <- df %>% 
            dplyr::group_by(Position) %>%
            dplyr::summarize(Mean=mean(Value), 
                                SD=sd(Value),
                                Upper=Mean+SD,
                                Lower=Mean-SD)
    # ADD GROUP ---
    dm$SampleID <- sampleid

    return(dm)
}




### GET FILES ---
files.matrix <- sapply(sampleids, function(x) { paste(dir.matrix, "/computematrix_", sprintf("%s", x), ".tsv.gz", sep="") })

### LOAD MATRIX DATA FOR EACH SAMPLE ---
list.dat <- list()
for(sampleid in sampleids){
    file.matrix <- files.matrix[sampleid]

    # LOAD DATA ---
    list.dat[[sampleid]] <- loadData(file.matrix, sampleid)

    cat("PROCESSED:", sampleid, "\n", sep=" ")
}


### GET MEAN PROFILE MATRIX FOR ALL SAMPLES ----
mat <- matrix(0, nrow=length(sampleids), ncol=60, dimnames=list(list.rds_atacseq_masterdata$metadata_mcrpc$Sample_Name, c(1:60)))

for(sampleid in sampleids){
    mat[sampleid,] <- as.numeric(list.dat[[sampleid]]$Mean)
}


### RESHAPE DATA ---
df <- reshape2::melt(as.matrix(mat))
colnames(df) <- c("SampleID","Position","Value")

### GET MEAN PROFILE ---
dm <- df %>% 
        dplyr::group_by(Position) %>%
        dplyr::summarize(Mean=mean(Value), 
                            SD=sd(Value),
                            Upper=Mean+SD,
                            Lower=Mean-SD)


### PLOT ---
p <- ggplot(data = dm, aes(x=Position, y=Mean)) +
            geom_ribbon(aes(ymin=Lower, ymax=Upper), fill="#d9d9d9", alpha=0.9) +        
            geom_line(aes(y=Mean), color="#E31A1C", alpha=0.8, size=0.6) +
            scale_x_continuous(breaks=c(1,11,21,30,40,50,60), labels=c("-3","-2","-1","TSS","1","2","3")) +
            scale_y_continuous(breaks=seq(0, 1.5, by=0.25), labels=seq(0, 1.5, by=0.25)) +
            coord_cartesian(xlim=c(1,60), ylim=c(0,1.5)) +
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
            ylab("ATAC-seq Profile") +
            xlab("Distance from TSS") + 
            ggtitle("")


# WRITE PLOT ---
file.plot <- file.path(dir.reproduce_fig, "supplementary_figure_03_f.pdf")
pdf(file.plot, height=2, width=2)
    grid.arrange(p, nrow=1, ncol=1) 
dev.off()
