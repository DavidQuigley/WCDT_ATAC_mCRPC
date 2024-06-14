###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### LOAD LIBRARIES ---
library("stringr")
library("data.table")
library("reshape2")
library("ggplot2")
library("gridExtra")

### DEFINE PATH ---
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/analysis/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

dir.tfpeaks2gene <- file.path(dir.reproduce_data, "tfpeaks2gene")
dir.footprints <- file.path(dir.reproduce_data, "tf_footprints/mcrpc_subtypes/data_footprints") 

### DEFINE FILES ---
file.tf2genes <- file.path(dir.tfpeaks2gene, "wcdt_tf2genes_r_0p4_pval_0p05_distTSS_500kb.tsv.gz")
file.tf_hits_rds <- file.path(dir.footprints, "tobias_mcrpc_subtypes_tf_hits.rds")
file.rds_genelist <- file.path(dir.reproduce_data, "genelist.rds")
file.utility_functions <- file.path(dir.reproduce_scripts, "00_prepare_data/00_01_utility_functions.R")

### LOAD FUNCTIONS ---
source(file.utility_functions)

### MOTIFS ---
motif_prefix_znf <- "MA0528.1_ZNF263"
motif_prefix_myc <- "MA0147.3_MYC"


###########################################################################################################################
### FUNCTION: get_tf_subtype() ---
get_tf_subtype <- function(dat_tfsubtypes, subtype){
    # EXTRACT TF BY SUBTYPE ---
    index_subtype <- which( colnames(dat_tfsubtypes) == subtype )
    index_tf <- which( dat_tfsubtypes[,index_subtype] >= 2 ) 
    dat_tfsubtypes <- dat_tfsubtypes[index_tf,]
    tfs_subtype <- dat_tfsubtypes$motif_prefix

    return(tfs_subtype)
}

### FUNCTION: getStats() ---
getStats <- function(list.feature_ids, genes_ref, tf){
    tf_query <- names(list.feature_ids)
    n_genes_ref <- length(genes_ref)
    n_genes_query <- as.numeric( unlist(lapply(list.feature_ids, function(x) length(x))) )
    n_ref_query <- as.numeric(unlist( lapply(list.feature_ids, function(x)  length( intersect(x, genes_ref) ) ) ) )
    percent_overlap <- round( (n_ref_query/n_genes_ref) * 100, 2)

    # PREPARE DATA ---
    dat <- data.frame(tf_ref=tf,
                    tf_query=tf_query,
                    n_genes_ref=n_genes_ref,
                    n_genes_query=n_genes_query,
                    n_ref_query=n_ref_query,
                    percent_overlap=percent_overlap)

    return(dat)
}




############################################################################################################################################
### LOAD TF HITS DATA ---
list.tf_hits <- readRDS(file=file.tf_hits_rds)
dat_tfsubtypes <- list.tf_hits$tf_footprints_hits_refined

### LOAD GENELIST ---
list.rds_genelist <- readRDS(file=file.rds_genelist)
genes_hallmark_myc <- unique( c( list.rds_genelist$msigdb_pathways$hallmark$HALLMARK_MYC_TARGETS_V1, list.rds_genelist$msigdb_pathways$hallmark$HALLMARK_MYC_TARGETS_V2) )

### LOAD TF2GENES DATA ---
dat_tf2genes <- data.table::fread(file=file.tf2genes, sep="\t", header=TRUE, nThread=50, data.table=FALSE, verbose=FALSE)
genes_myc <- unique(dat_tf2genes$Gene[which(dat_tf2genes$motif_prefix == motif_prefix_myc)])



############################################################################################################################################
### GET TFs ASSOCIATED WITH DNPC (AR-NE-) ---
tfs_dnpc <- get_tf_subtype(dat_tfsubtypes, subtype="ARnNEn")


### EXTRACT TF-GENES DATA ---
list.tf_targets <- list()
for(i in 1:length(tfs_dnpc)){
    motif_prefix <- tfs_dnpc[i]
    tf <- stringr::str_split(tfs_dnpc[i], "_")[[1]][2]
    list.tf_targets[[tf]] <- dat_tf2genes[which(dat_tf2genes$motif_prefix == motif_prefix),]
}

### GET TARGET GENES ---
list.feature_ids <- lapply(list.tf_targets, function(x) unique(x$Gene) )

### GET GENES OVERLAP STATS ---
df_stats <- getStats(list.feature_ids, genes_ref=genes_myc, tf="MYC")
df_stats_hallmark <- getStats(list.feature_ids, genes_ref=genes_hallmark_myc, tf="MYC_sig")

### MERGE DATA ---
df_stats$percent_overlap_hallmark <- df_stats_hallmark$percent_overlap

### ORDER DATA BY PERCENT OVERLAP ---
df_stats <- df_stats[order(df_stats$percent_overlap, decreasing=TRUE),]



############################################################################################################################################
### ADD GROUP INFO --
df_stats$group <- ifelse(df_stats$tf_query == "ZNF263", "HIGHLIGHT", "BG")

### FACTORIZE DATA ---
df_stats$tf_query <- factor(df_stats$tf_query, levels=df_stats$tf_query)
df_stats$group <- factor(df_stats$group, levels=c("HIGHLIGHT","BG"))



### DEFINE AXES: PRIMARY & SECONDARY ----------------
ylim.prim <- c(40, 100)
ylim.sec <- c(20, 80)

b <- diff(ylim.prim)/diff(ylim.sec)
a <- ylim.prim[1] - (b*ylim.sec[1])


### PLOT BAR-LINE DOUBLE AXIS ----
p1 <- ggplot(df_stats, aes(x=tf_query, y=percent_overlap)) +
            geom_bar(aes(fill=group), stat="identity", color="#000000", width=0.8, size=0.2) +
            #geom_line(stat="identity", color="red", size=0.5, alpha=0.7) +
            geom_line(aes(y = a+ (percent_overlap_hallmark*b)), color="#e41a1c", size=0.5, alpha=0.7, group = 1) +
            geom_point(aes(y = a+ (percent_overlap_hallmark*b)), fill="#e41a1c", color="#e41a1c", shape=21, stroke=0.25, size=0.7, alpha=0.7, na.rm=TRUE) +
            coord_cartesian(ylim=c(40,100)) +
            scale_y_continuous(breaks=seq(40, 100, by=20), sec.axis = sec_axis(~ (. - a)/b, name = "% of target genes overlap\n with Hallmark MYC Sig. pathway genes") ) +
            scale_fill_manual(values=c("#ff7f00","#006094")) +
            theme(
               axis.text.x = element_text(size = 4, color="#000000", angle=90, hjust=1, vjust=0.5),
               axis.text.y = element_text(size = 4, color="#000000"),
               axis.text.y.right = element_text(size = 4, color="#e41a1c"),
               axis.title = element_text(size = 4, color="#000000"),
               axis.title.y.right = element_text(size = 4, color="#e41a1c"),
               plot.title = element_text(size = 4, color="#000000", hjust=0.5),
               panel.grid.major.y = element_blank(),
               panel.grid.major.x = element_blank(),
               panel.grid.minor = element_blank(),
               axis.ticks.y = element_line(size=0.2, color="#000000"),
               axis.ticks.x = element_line(size=0.2, color="#000000"),
               axis.line.y.right = element_line(color="#e41a1c"), 
               axis.ticks.y.right = element_line(color="#e41a1c"),               
               strip.text = element_text(size=4, color="#000000"),
               strip.background = element_rect(fill="#FFFFFF", color="#FFFFFF"),
               panel.background = element_rect(fill="#FFFFFF", color="#000000"),
               legend.text = element_text(size = 4, color="#000000"),
               legend.title = element_blank(),
               legend.key.size = unit(0.5, "cm"),
               legend.position = "none") +
            xlab("") +
            ylab("% of target genes overlap\n with that of MYC") + 
            ggtitle("")


### GENERATE PLOT ---
file.plot <- file.path(dir.reproduce_fig, "figure_6_c.pdf")
pdf(file.plot, height=1.5, width=2.5)
    grid.arrange(p1, nrow=1, ncol=1)
dev.off()

