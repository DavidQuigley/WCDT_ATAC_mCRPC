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
library("ggrepel")

### DEFINE PATH ---
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

dir.footprints <- file.path(dir.reproduce_data, "tf_footprints/mcrpc_subtypes")


### DEFINE FILES ---
file.motifdb <- file.path(dir.reproduce, "reference/JASPAR2018_motif_gene_conversion.tsv")
file.rds_bindetect_results <- file.path(dir.footprints, "data_footprints/tobias_mcrpc_subtypes_bindetect_results.rds")

### SUBTYPES ---
subtypes <- c("ARpNEn","ARpNEp","ARlNEn","ARnNEp","ARnNEn")
cpalette <- c("#e31a1c","#ff7f00","#fb9a99","#33a02c","#1f78b4")

### GET PAIRWISE ALL GROUP COMBINATIONS ---
dm_grp <- do.call(rbind.data.frame, combn(x=subtypes, m=2, simplify=FALSE)) %>%
                dplyr::mutate_all(as.character)
colnames(dm_grp) <- c("Group1","Group2")
rownames(dm_grp) <- NULL
dm_grp$id <- apply(dm_grp, 1, function(x) paste(x[1], x[2], sep="_"))
dm_grp$f_name <- apply(dm_grp, 1, function(x) paste("tobias_diff", paste(x[1], x[2], sep="_"), sep="_"))

### ADD COLOR ---
dm_grp$Color1 <- NA
dm_grp$Color2 <- NA
for(i in 1:length(cpalette)){
    dm_grp$Color1[which(dm_grp$Group1 == subtypes[i])] <- cpalette[i]
    dm_grp$Color2[which(dm_grp$Group2 == subtypes[i])] <- cpalette[i]
}


### FUNCTION: ParseMotifDB() ---
ParseMotifDB <- function(file.motifdb){
    # LOAD MOTIFDB DATA ---
    motifdb <- read.delim(file.motifdb, header=TRUE, stringsAsFactors=FALSE)
    motifdb <- subset(motifdb, motifdb$species == "Homo sapiens")
    
    motifdb$motif_prefix <- stringr::str_replace_all(motifdb$motif_prefix, "::", "")
    motifdb$motif_prefix <- stringr::str_replace_all(motifdb$motif_prefix, "[(]", "")
    motifdb$motif_prefix <- stringr::str_replace_all(motifdb$motif_prefix, "[)]", "")

    motif_prefix <- stringr::str_replace_all(motifdb$motif_name, "::", "_")
    motif_prefix <- stringr::str_replace_all(motif_prefix, "[(]", "_")
    motif_prefix <- stringr::str_replace_all(motif_prefix, "[)]", "")
    motif_prefix <- stringr::str_replace_all(motif_prefix, "[.]", "_")

    motifdb$motif_prefix_2 <- paste(motifdb$motif_id, motif_prefix, sep="_")
    motifdb <- motifdb[order(motifdb$motif_id, decreasing=FALSE),]

    rownames(motifdb) <- NULL

    motifdb <- subset(motifdb, select=c("motif_id","motif_name","motif_prefix_2"))
    colnames(motifdb) <- c("motif_id","motif_name","motif_prefix")
    motifdb <- motifdb[!duplicated(motifdb),]

    return(motifdb)
}


### FUNCTION: parseDataSummary() ---
parseDataSummary <- function(dat, motif_db){
    # ADD MOTIF PREFIX ---
    for(i in 1:nrow(dat)){
        dat$motif_prefix[i] <- motif_db$motif_prefix[which(motif_db$motif_id == dat$motif_id[i])]
        dat$name[i] <- motif_db$motif_name[which(motif_db$motif_id == dat$motif_id[i])]
    }

    # GET COLUMNS OF INTEREST ---
    y_change <- colnames(dat)[stringr::str_detect(colnames(dat), "_change")]
    y_pvalue <- colnames(dat)[stringr::str_detect(colnames(dat), "_pvalue")]

    # PREPARE DATA ---
    df <- subset(dat, select=c("motif_prefix","name","motif_id",y_change,y_pvalue))
    colnames(df) <- c("motif_prefix","motif_name","motif_id","score","pvalue")

    # CORRECT PVALUE ---
    index0 <- which(df$pvalue == 0)
    if(length(index0) != 0){
        df$pvalue[index0] <- .Machine$double.eps #smallest value
    }

    # PVALUE CORRECTION ---
    df$fdr <- p.adjust(df$pvalue, method="fdr", n = length(df$pvalue))
    df$nlogfdr <- -log10(df$fdr)

    return(df)
}


### FUNCTION:getParams()  ---
getParams <- function(df, id){
    y_cutoff <- as.numeric(round(quantile(df$nlogfdr, 0.8), 0))
    xp_cutoff <- as.numeric(round(quantile( df$score , 0.85), 1))
    xn_cutoff <- as.numeric(round(quantile( df$score , 0.15), 1))
    
    score_threshold <- switch(id,  
                        "ARpNEn_ARpNEp" = 0.05,
                        "ARpNEn_ARlNEn" = 0.1,
                        "ARpNEn_ARnNEp" = 0.1,
                        "ARpNEn_ARnNEn" = 0.1,
                        "ARpNEp_ARlNEn" = 0.1,
                        "ARpNEp_ARnNEp" = 0.1,
                        "ARpNEp_ARnNEn" = 0.1,
                        "ARlNEn_ARnNEp" = 0.1,
                        "ARlNEn_ARnNEn" = 0.1,
                        "ARnNEp_ARnNEn" = 0.1)

    score_threshold_p <- switch(id,  
                        "ARpNEn_ARpNEp" = xp_cutoff,
                        "ARpNEn_ARlNEn" = xp_cutoff,
                        "ARpNEn_ARnNEp" = xp_cutoff,
                        "ARpNEn_ARnNEn" = xp_cutoff,
                        "ARpNEp_ARlNEn" = xp_cutoff,
                        "ARpNEp_ARnNEp" = xp_cutoff,
                        "ARpNEp_ARnNEn" = xp_cutoff,
                        "ARlNEn_ARnNEp" = xp_cutoff,
                        "ARlNEn_ARnNEn" = xp_cutoff,
                        "ARnNEp_ARnNEn" = xp_cutoff)

    score_threshold_n <- switch(id,  
                        "ARpNEn_ARpNEp" = xn_cutoff,
                        "ARpNEn_ARlNEn" = xn_cutoff,
                        "ARpNEn_ARnNEp" = xn_cutoff,
                        "ARpNEn_ARnNEn" = xn_cutoff,
                        "ARpNEp_ARlNEn" = xn_cutoff,
                        "ARpNEp_ARnNEp" = xn_cutoff,
                        "ARpNEp_ARnNEn" = xn_cutoff,
                        "ARlNEn_ARnNEp" = xn_cutoff,
                        "ARlNEn_ARnNEn" = xn_cutoff,
                        "ARnNEp_ARnNEn" = xn_cutoff)

    score_upperbound <- switch(id,  
                        "ARpNEn_ARpNEp" = 0.4,
                        "ARpNEn_ARlNEn" = 0.6,
                        "ARpNEn_ARnNEp" = 0.6,
                        "ARpNEn_ARnNEn" = 0.6,
                        "ARpNEp_ARlNEn" = 0.4,
                        "ARpNEp_ARnNEp" = 0.4,
                        "ARpNEp_ARnNEn" = 0.6,
                        "ARlNEn_ARnNEp" = 0.6,
                        "ARlNEn_ARnNEn" = 0.6,
                        "ARnNEp_ARnNEn" = 0.6)

    nlogp_threshold <- switch(id,  
                        "ARpNEn_ARpNEp" = y_cutoff,
                        "ARpNEn_ARlNEn" = y_cutoff,
                        "ARpNEn_ARnNEp" = y_cutoff,
                        "ARpNEn_ARnNEn" = y_cutoff,
                        "ARpNEp_ARlNEn" = y_cutoff,
                        "ARpNEp_ARnNEp" = y_cutoff,
                        "ARpNEp_ARnNEn" = y_cutoff,
                        "ARlNEn_ARnNEp" = y_cutoff,
                        "ARlNEn_ARnNEn" = y_cutoff,
                        "ARnNEp_ARnNEn" = y_cutoff)

    nlogp_upperbound <- switch(id,  
                        "ARpNEn_ARpNEp" = 180,
                        "ARpNEn_ARlNEn" = 220,
                        "ARpNEn_ARnNEp" = 200,
                        "ARpNEn_ARnNEn" = 240,
                        "ARpNEp_ARlNEn" = 220,
                        "ARpNEp_ARnNEp" = 200,
                        "ARpNEp_ARnNEn" = 220,
                        "ARlNEn_ARnNEp" = 200,
                        "ARlNEn_ARnNEn" = 220,
                        "ARnNEp_ARnNEn" = 220)

    list.output <- list(score_threshold_pos=score_threshold,
                        score_threshold_neg=-score_threshold,
                        score_upperbound=score_upperbound,
                        nlogp_threshold=nlogp_threshold,
                        nlogp_upperbound=nlogp_upperbound)

    return(list.output)
}




### FUNCTION: getPlot() ---
getPlot <- function(df, score.threshold_pos, score.threshold_neg, nlogp.threshold, score.upperbound, nlogp.upperbound, colors_pack, compr.name, labels=FALSE){
    yval <- nlogp.threshold
    xval_pos <- score.threshold_pos
    xval_neg <- score.threshold_neg
    cbPalette <- c(colors_pack, "#98a3a5")

    # GROUP TF-MOTIF ---
    df$Group <- "GROUP_3"
    df$Group[which((df$score >= score.threshold_pos) & (df$nlogfdr >= nlogp.threshold))] <- "GROUP_1"
    df$Group[which((df$score <= score.threshold_neg) & (df$nlogfdr >= nlogp.threshold))] <- "GROUP_2"

    # ADD TF MOTIF NAMES TO LABEL ---
    df$Motif_Label <- NA
    if(labels){
        df$Motif_Label[which(df$Group == "GROUP_1")] <- df$motif_name[which(df$Group == "GROUP_1")]
        df$Motif_Label[which(df$Group == "GROUP_2")] <- df$motif_name[which(df$Group == "GROUP_2")]
    }

    # ADJUST FOLDCHANGE FOR VIZ ---
    if(max(df$score) > score.upperbound){
        df$score[which(df$score > score.upperbound)] <- score.upperbound 
    }

    if(min(df$score) < -score.upperbound){
        df$score[which(df$score < -score.upperbound)] <- -score.upperbound
    }

    # ADJUST nlogfdr FOR VIZ ---
    if(max(df$nlogfdr) > nlogp.upperbound){
        df$nlogfdr[which(df$nlogfdr > nlogp.upperbound)] <- nlogp.upperbound 
    }

    # ORDER DATA ---
    df <- df[order(df$Group, decreasing=TRUE),]

    # FACTORIZE ---
    df$Group <- factor(df$Group, levels=c("GROUP_1","GROUP_2","GROUP_3"))

    # PLOT ---
	p <- ggplot(df, aes(x=score, y=nlogfdr, label=Motif_Label)) +         
			geom_hline(yintercept = yval, colour="gray70", linetype=4, alpha=0.7, size=0.25) +
			geom_vline(xintercept = xval_pos, colour="gray70", linetype=4, alpha=0.7, size=0.25) +		
			geom_vline(xintercept = xval_neg, colour="gray70", linetype=4, alpha=0.7, size=0.25) +
			#geom_point(aes(fill=Group, color=Group), shape = 21, stroke=0.3, size=0.8, alpha=0.8) +
            geom_point(aes(fill=Group, color=Group, size=Group), shape = 21, stroke=0.3, alpha=0.8) +
			scale_fill_manual(values=cbPalette) +
			scale_color_manual(values=cbPalette) +
            scale_size_manual(values=c(0.8,0.8,0.2)) +
            scale_x_continuous(breaks=round(seq(-score.upperbound, score.upperbound, 0.2), 1), 
                                labels=round(seq(-score.upperbound, score.upperbound, 0.2), 1) ) +
            scale_y_continuous(breaks=seq(0, nlogp.upperbound, 50)) +
			coord_cartesian(xlim=c(-score.upperbound, score.upperbound), ylim=c(0, nlogp.upperbound)) +	
			theme(
				axis.text = element_text(size = 5, color="#000000"),
				axis.title = element_text(size = 5, color="#000000"),
				strip.text = element_text(size = 5, color="#000000"),
				strip.background = element_rect(fill = "#FFFFFF", colour = "#FFFFFF"),
				plot.title = element_text(size = 5, color="#000000", hjust=0.5),
				panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				axis.ticks = element_line(size=0.2, color="#000000"),	
				panel.background = element_rect(fill="#FFFFFF", color="#000000"),
				legend.text = element_text(size = 5, color="#000000"),
			    legend.title = element_blank(),
			    legend.key.size = unit(0.3, "cm"),			
			    legend.position="none") +
			xlab("Differential Binding Score") + 
			ylab("-log10(FDR pvalue)") + 
			ggtitle(compr.name)     

    # ADD GEOME TEXT ---
    if(labels){
        #p <- p + geom_text(aes(label=Motif_Label, alpha=0.3, hjust=1, vjust=-1, lineheight=1.5), color="#000000", size=1.5)
        p <- p + ggrepel::geom_text_repel(aes(label=Motif_Label), segment.color=NA, color="#000000", alpha=0.8, size=1.5, max.overlaps=25)
    }

    return(p)
}

### GET MOTIFS ---
motif_db <- ParseMotifDB(file.motifdb)

### LOAD BINDETECT RESULTS ---
list.bindetect_results <- readRDS(file=file.rds_bindetect_results)




### LOOP FOR EACH PAIR-WISE COMPARISON ---
list.plotVolcano <- list()
for(k in 1:nrow(dm_grp)){
    id <- dm_grp$id[k]
    f_name <- dm_grp$f_name[k]
    group_id_1 <- dm_grp$Group1[k]
    group_id_2 <- dm_grp$Group2[k]
    color_1 <- dm_grp$Color1[k]
    color_2 <- dm_grp$Color2[k]
    colors_pack <- c(color_1, color_2)

    cat("START:", id, "\n", sep=" ")  

    # GET BINDETECT REULTS DATA ---
    dat <- list.bindetect_results[[id]]

    # PARSE DATA ---
    df <- parseDataSummary(dat, motif_db)

    # GET PARAMS ---
    list.params <- getParams(df, id)
    score_threshold_pos <- list.params$score_threshold_pos
    score_threshold_neg <- list.params$score_threshold_neg
    score_upperbound <- list.params$score_upperbound
    nlogp_threshold <- list.params$nlogp_threshold
    nlogp_upperbound <- list.params$nlogp_upperbound

    # GET PLOT ---
    list.plotVolcano[[id]] <- getPlot(df=df, score.threshold_pos=score_threshold_pos, score.threshold_neg=score_threshold_neg, nlogp.threshold=nlogp_threshold, score.upperbound=score_upperbound, nlogp.upperbound=nlogp_upperbound, colors_pack, compr.name=id, labels=FALSE)

    cat("PROCESSED:", id, "\n", sep=" ")  
}

### GENERATE PLOT ---
file.plot <- file.path(dir.reproduce_fig, "supplementary_figure_15.pdf")
pdf(file.plot, height=1.8, width=1.8)
for(k in 1:length(list.plotVolcano)){    
    gridExtra::grid.arrange(list.plotVolcano[[k]], nrow=1, ncol=1)
}
dev.off()
