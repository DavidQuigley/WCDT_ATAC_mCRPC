###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### LOAD LIBRARIES ---
library("stringr")
library("data.table")
library("GenomicRanges")
library("IRanges")
library("ggplot2")
library("gridExtra")
library("Cairo")


### DEFINE PATH ---
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/analysis/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

dir.tfpeaks2gene <- file.path(dir.reproduce_data, "tfpeaks2gene")

### DEFINE FILES ---
file.rds_chipseq_bed <- file.path(dir.reproduce_data, "chipseq_bed.rds")

list.file_dat <- list()
list.file_dat[[1]] <- file.path(dir.tfpeaks2gene, "wcdt_tfpeaks2genes_correlated_r_0p4_pval_0p05.tsv.gz")
list.file_dat[[2]] <- file.path(dir.tfpeaks2gene, "wcdt_tfpeaks2genes_correlated_r_0p5_pval_0p05.tsv.gz")
list.file_dat[[3]] <- file.path(dir.tfpeaks2gene, "wcdt_tfpeaks2genes_correlated_r_0p5_pval_0p01.tsv.gz")
list.file_dat[[4]] <- file.path(dir.tfpeaks2gene, "wcdt_tfpeaks2genes_correlated_r_0p5_fdr_0p05.tsv.gz")
list.file_dat[[5]] <- file.path(dir.tfpeaks2gene, "wcdt_tfpeaks2genes_correlated_r_0p5_fdr_0p01.tsv.gz")
list.file_dat[[6]] <- file.path(dir.tfpeaks2gene, "wcdt_tfpeaks2genes_correlated_r_0p6_pval_0p05.tsv.gz")
list.file_dat[[7]] <- file.path(dir.tfpeaks2gene, "wcdt_tfpeaks2genes_correlated_r_0p7_pval_0p05.tsv.gz")
list.file_dat[[8]] <- file.path(dir.tfpeaks2gene, "wcdt_tfpeaks2genes_correlated_r_0p6_fdr_0p05.tsv.gz")
list.file_dat[[9]] <- file.path(dir.tfpeaks2gene, "wcdt_tfpeaks2genes_correlated_r_0p7_fdr_0p05.tsv.gz")
list.file_dat[[10]] <- file.path(dir.tfpeaks2gene, "wcdt_tfpeaks2genes_correlated_r_0p4_pval_0p05_distTSS_500kb.tsv.gz")
list.file_dat[[11]] <- file.path(dir.tfpeaks2gene, "wcdt_tfpeaks2genes_correlated_r_0p4_pval_0p05_distTSS_1000kb.tsv.gz")



### FUNCTION: loadData() ---
loadData <- function(file.dat){
    # LOAD TFPEAKS2GENE DATA ---
    dat <- data.table::fread(file=file.dat, sep="\t", header=TRUE, nThread=60, data.table=FALSE, verbose=FALSE)
    dat$chr <- unlist(lapply(stringr::str_split(dat$Peak, "_"), function(x) x[1]))
    dat$start <- as.numeric(unlist(lapply(stringr::str_split(dat$Peak, "_"), function(x) x[2])))
    dat$end <- as.numeric(unlist(lapply(stringr::str_split(dat$Peak, "_"), function(x) x[3])))
    gr <- GenomicRanges::makeGRangesFromDataFrame(df=dat, keep.extra.columns=TRUE)
    return(gr)
}

### FUNCTION: compute_overlaps() ---
compute_overlaps <- function(gr, gr_overlap, group_label){
    # ADD H3K27ac OVERLAP TAG ---
    y <- rep(0, length(gr))
    y[which(gr$Peak %in% gr_overlap$Peak)] <- 1
    gr$H3K27ac <- y
    gr$Group <- group_label

    # ADD GROUP INFO ---
    dm <- gr %>%
            as.data.frame() %>%
            dplyr::count(Group, H3K27ac, sort=FALSE, name="Freq") %>%
            dplyr::mutate(Percent=round( (Freq/sum(Freq)) * 100 , 2)  )
    
    return(dm)
}

### FUNCTION: get_overlaps() ---
get_overlaps <- function(gr_a, gr_anb){
    y_anb <- length(gr_anb)
    y_a <- length(gr_a)
    p <- round((y_anb/y_a) * 100, 2)
    return(p)
}


################################################################################################################
### LOAD CHIPSEQ DATA ---
list.rds_chipseq_bed <- readRDS(file=file.rds_chipseq_bed)
gr_h3k27ac <- list.rds_chipseq_bed$chipseq_h3k27ac_kron



### GET LOAD DATA AND CONVERT TO GRANGES ---
list.gr <- list()
for(i in 1:length(list.file_dat)){
    list.gr[[i]] <- loadData(file.dat=list.file_dat[[i]])
}


### GET OVERLAP ---
list.gr_overlap <- list()
for(i in 1:length(list.file_dat)){
    list.gr_overlap[[i]] <- IRanges::subsetByOverlaps(x=list.gr[[i]], ranges=gr_h3k27ac, type="any")
}




### COMPUTE H3K27ac OVERLAP PERCENTAGE ---
list.df <- list()
for(i in c(1:length(list.file_dat))){
    percent <- get_overlaps(gr_a=list.gr[[i]], gr_anb=list.gr_overlap[[i]])
    list.df[[i]] <- data.frame(Group=LETTERS[i], Percent=percent)
}

### AGGREGATE DATA ---
df <- do.call(rbind.data.frame, list.df)
rownames(df) <- NULL

### ADD TAG ---
df$Tag <- ifelse(df$Group == "J", "PASS", "FAIL")

### FACTORIZE ---
df$Group <- factor(df$Group, levels=LETTERS[1:11])
df$Tag <- factor(df$Tag, levels=c("PASS","FAIL"))

### UNICODE
# <= \u2264
# >= \u2265



### PLOT: PERCENT ---
p <- ggplot(df, aes(x=Group, y=Percent)) + 
        geom_bar(aes(fill=Tag), stat="identity", color="#000000", width=0.8, alpha=0.7, size=0.5) + 
        scale_fill_manual(values=c("#ff7f00","#006094")) +
        coord_cartesian(ylim=c(60,75)) +
        scale_y_continuous(breaks=seq(60,75,by=5)) +
        scale_x_discrete(labels=c("A" = "R \u2265 |0.4|, p \u2264 0.05",
                                "B" = "R \u2265 |0.5|, p \u2264 0.05",
                                "C" = "R \u2265 |0.5|, p \u2264 0.01",
                                "D" = "R \u2265 |0.5|, q \u2264 0.05",
                                "E" = "R \u2265 |0.5|, q \u2264 0.01",
                                "F" = "R \u2265 |0.6|, p \u2264 0.05",
                                "G" = "R \u2265 |0.7|, p \u2264 0.05",
                                "H" = "R \u2265 |0.6|, q \u2264 0.05",
                                "I" = "R \u2265 |0.7|, q \u2264 0.05",
                                "J" = "R \u2265 |0.4|, p \u2264 0.05, d \u2264 500kb",
                                "K" = "R \u2265 |0.4|, p \u2264 0.05, d \u2264 1000kb")) +
        theme(
          axis.text.x = element_text(size = 5, color="#000000", angle=45, hjust=1, vjust=1),
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
          legend.position = "none") +
      ylab("Percent of ATAC-seq peaks \n overlapped with H3K27ac marks") +
      xlab("") + 
      ggtitle("") 


### GENERATE PLOT ---
file.plot <- file.path(dir.reproduce_fig, "supplementary_figure_20_a.pdf")
Cairo::CairoPDF(file=file.plot, height=2.25, width=2)
    grid.arrange(p, nrow=1, ncol=1)
dev.off()

