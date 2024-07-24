###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### DEFINE LIBRARIES ---
library("stringr")
library("data.table")
library("dplyr")
library("reshape2")
library("ggplot2")
library("gridExtra")

### DEFINE PATH ---
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

dir.base <- file.path("/data1/datasets_1/human_prostate_WCDT/atac_solid")
dir.process <- file.path(dir.base, "processed/2020_09_22/processed")
dir.output <- tempdir() # REPLACE PATH
dir.scripts <- tempdir() # REPLACE PATH
dir.logs <- tempdir() # REPLACE PATH


### DEFINE FILES ----
file.rds_atacseq_masterdata <- file.path(dir.reproduce_data, "wcdt_atacseq_masterdata.rds")
file.script <- file.path(dir.reproduce_scripts, "function_compute_insert_sizes.R")


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

### LOAD ATAC-SEQ MASTER DATA ---
list.rds_atacseq_masterdata <- readRDS(file=file.rds_atacseq_masterdata)
sampleids <- list.rds_atacseq_masterdata$metadata_mcrpc$Sample_Name

### GET BAM FILES ---
files.bam <- sapply(sampleids, function(x) { paste(dir.process, "/", sprintf("%s/alignment/%s", x, x), ".bowtie.sorted.nodup.shifted.bam", sep="") })

### LOOP FOR EACH SAMPLE ---
list.sh <- list()
for(sampleid in sampleids){
    file.bam <- files.bam[sampleid]
    file.output <- file.path(dir.output, paste("insertsizes_", sampleid, ".tsv", sep="" ))

    # GET COMMANDS ---
    list.cmd_sge <- get.SGEcmd(jobid=paste("fs", sampleid, sep="_"), dir.logs)
    cmd_r <- paste("Rscript", file.script, 
                    "-b", file.bam,
                    "-o", file.output,
                    "-s", sampleid, sep=" ")

    # MERGE COMMANDS ---
    cmd <- c(unlist(list.cmd_sge), cmd_r)

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

### GET FILES ---
files.dat <- sapply(sampleids, function(x) { paste(dir.output, "/insertsizes_", sprintf("%s", x), ".tsv", sep="") })

### LOAD DATA FOR EACH SAMPLE ---
list.dat <- list()
for(sampleid in sampleids){
    file.dat <- files.dat[sampleid]

    # LOAD DATA ---
    dat <- data.table::fread(file=file.dat, sep="\t", header=TRUE, nThread=55,  data.table=FALSE, verbose=FALSE)

    # LOAD DATA ---
    list.dat[[sampleid]] <- dat
}


### AGGREAGTE DATA ---
df <- do.call(rbind.data.frame, list.dat)
rownames(df) <- NULL

df$Count <- NULL
df$CountLog <- NULL

### CONVERT TO MATRIX ---
dt_df <- data.table::setDT(x=df)
dt_mat <- as.data.frame(data.table::dcast(dt_df, SampleID ~ InsertSizes, value.var = "Proportion"))
dt_mat[is.na(dt_mat)] <- 0


### RESHAPE DATA ---
df_mat <- reshape2::melt(dt_mat, id.vars="SampleID", variable.name="InsertSize", value.name="Percent")

### GET MEAN PROFILE ---
dm <- df_mat %>% 
        dplyr::group_by(InsertSize) %>%
        dplyr::summarize(Mean=mean(Percent), 
                            SD=sd(Percent),
                            Upper=abs(Mean+SD),
                            Lower=abs(Mean-SD))

dm$InsertSize <- as.numeric(as.character(dm$InsertSize))



### GENERATE PLOT ---
p <- ggplot(dm, aes(x=InsertSize, y=Mean)) +
            geom_ribbon(aes(ymin=Lower, ymax=Upper), fill="#d9d9d9", alpha=0.9) +        
            geom_line(aes(y=Mean), color="#E31A1C", alpha=0.8, size=0.6) +
            scale_x_continuous(breaks=seq(0, 500, by=100), labels=seq(0, 500, by=100)) +
            scale_y_continuous(breaks=seq(0, 1.5, by=0.25), labels=seq(0, 1.5, by=0.25)) +
            coord_cartesian(xlim=c(0,500), ylim=c(0,1.5)) +
            theme(
                axis.text.x = element_text(size = 6, color="#000000"),
                axis.text.y = element_text(size = 6, color="#000000"),
                axis.title = element_text(size = 6, color="#000000"),
                plot.title = element_text(size = 6, color="#000000", hjust=0.5),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.ticks = element_line(size=0.2, color="#000000"), 
                strip.text = element_text(size=6, color="#000000"),
                strip.background = element_rect(fill="#FFFFFF", color="#FFFFFF"),
                panel.background = element_rect(fill="#FFFFFF", color="#000000"),
                legend.text = element_text(size = 6, color="#000000"),
                legend.title = element_blank(),
                legend.key.size = unit(0.3, "cm"),
                legend.position = "none") +
            ylab("Fragments (%)") +
            xlab("Fragment Length (bp)") + 
            ggtitle("")

# WRITE PLOT ---
file.plot <- file.path(dir.reproduce_fig, "supplementary_figure_03_g.pdf")
pdf(file.plot, height=2, width=2)
    grid.arrange(p, nrow=1, ncol=1) 
dev.off()
