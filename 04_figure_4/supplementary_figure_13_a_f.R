###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### LOAD LIBRARIES ---
library("stringr")
library("data.table")
library("pracma")
library("dplyr")
library("reshape2")
library("ggplot2")
library("gridExtra")

### DEFINE PATH ---
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/analysis/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

dir.footprints <- file.path(dir.reproduce_data, "tf_footprints/mcrpc_subtypes")
dir.bam <- file.path(dir.reproduce_data, "tf_footprints/mcrpc_subtypes/bam")
dir.bed <- file.path(dir.reproduce_data, "tf_footprints/mcrpc_subtypes/bed")

### DEFINE FILES ---
file.rds_footprints <- file.path(dir.footprints, "data_footprints/expt_bound_pvalue_footprint_summary.rds")
file.performance_metrics <- file.path(dir.footprints, "data_footprints/expt_footprint_chipseq_prediction_performance_metrics.tsv")


### SUBTYPES ---
subtypes <- c("ARpNEn","ARlNEn","ARpNEp","ARnNEp","ARnNEn")

### PVALUES ---
pvals <- c(0.5, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001)
names(pvals) <- c("p5","p1","p05","p01","p005","p001","p0005","p0001")

### MOTIFS ---
motifs <- c("ExtendedSite_AR","MA0148.3_FOXA1","MA0901.1_HOXB13")
color_palette <- c("#e31a1c","#fb9a99","#ff7f00","#33a02c","#1f78b4")

### FUNCTION: prepareData() ---
prepareData <- function(df,  motif, subtype){
    # CREATE LOWER BOUND AND UPPER BOUND ---
    df0 <- data.frame(motif=motif, subtype=subtype, pval="p0", sensitivity=0, fpr=0, ppv=0)
    df1 <- data.frame(motif=motif, subtype=subtype, pval="p11", sensitivity=1, fpr=1, ppv=1)

    # COMPILE DATA ---
    dm <- rbind(df0, df, df1)

    return(dm)
}


### FUNCTION: getROC() ---
getROC <- function(dm, motif, color_palette, text_label=FALSE){
    title_name <- motif

    # PLOT ---
    p1 <- ggplot(dm, aes(x=fpr, y=sensitivity, label=pval)) +
            geom_line(aes(color=subtype), size=0.5, alpha=0.7) +
            geom_point(aes(color=subtype, fill=subtype), shape=21, stroke=0.25, size=1, alpha=0.7, na.rm=TRUE) +
            geom_abline(intercept=0, linetype=2, color="#969696") +	
            coord_cartesian(xlim=c(0,1), ylim=c(0,1) ) +
            scale_x_continuous(breaks=seq(0, 1, by=0.2) ) +
            scale_y_continuous(breaks=seq(0, 1, by=0.2) ) +
            scale_fill_manual(values=color_palette) +
            scale_color_manual(values=color_palette) +
            theme(
               axis.text.x = element_text(size = 6, color="#000000"),
               axis.text.y = element_text(size = 6, color="#000000"),
               axis.title = element_text(size = 6, color="#000000"),
               plot.title = element_text(size = 8, color="#000000", hjust=0.5),
               panel.grid.major.y = element_blank(),
               panel.grid.major.x = element_blank(),
               panel.grid.minor = element_blank(),
               axis.ticks.y = element_line(size=0.2, color="#000000"),
               axis.ticks.x = element_line(size=0.2, color="#000000"),
               strip.text = element_text(size=7, color="#000000"),
               strip.background = element_rect(fill="#FFFFFF", color="#FFFFFF"),
               panel.background = element_rect(fill="#FFFFFF", color="#000000"),
               legend.text = element_text(size = 5, color="#000000"),
               legend.title = element_blank(),
               legend.key.size = unit(0.4, "cm"),
               legend.position = "none") +
            xlab("False Positive Rate") +
            ylab("True Positive Rate (Sensitivity)") + 
            ggtitle(title_name)

    if(text_label){
        p1 <- p1 + ggrepel::geom_label_repel(aes(label=pval), segment.color="#BDBDBD", color="#000000", alpha=0.8, size=1.5, max.overlaps=Inf, na.rm=TRUE)
    }

    return(p1)            
}

### FUNCTION: getAUC() ---
getAUC <- function(dm, subtypes, motif){
    list.res <- list()
    for(subtype in subtypes){
        d <- dm[which(dm$subtype == subtype),]
        d <- d[order(d$fpr, decreasing=FALSE),]

        # COMPUTE AUC ---
        auc <- pracma::trapz(x=d$fpr, y=d$sensitivity)

        # STORE RESULT ---
        list.res[[subtype]] <- data.frame(motif=motif, subtype=subtype, auc=auc)
    }

    # AGGREGATE RESULT ---
    df <- do.call(rbind.data.frame, list.res)
    rownames(df) <- NULL

    return(df)
}

### FUNCTION: plotBar() ----
plotBar <- function(df, subtypes, color_palette, motif){
    # FACTORIZE ---
    df$subtype <- factor(df$subtype, levels=subtypes)

    # PLOT ---
    p <- ggplot(df, aes(x=subtype, y=auc)) +
            geom_bar(aes(fill=subtype), stat="identity") +
            coord_cartesian(ylim=c(0,0.8)) +
            scale_y_continuous(breaks=seq(0,0.8,by=0.1)) +
            scale_fill_manual(values=color_palette) +
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
                legend.position = "none") +
            xlab("") +
            ylab("AUC") + 
            ggtitle(motif) 
    
    return(p)
}


##################################################################################################
### LOAD PERFORMANCE METRIC ---
dat <- data.table::fread(file=file.performance_metrics, sep="\t", header=TRUE, nThread=1, data.table=FALSE, verbose=FALSE)



### PREPARE DATA ---
list.dm <- list()
for(motif in motifs){
    # SUBSET BY MOTIF ---
    dat1 <- dat[which(dat$motif == motif),]

    # LOOP BY SUBTYPE ---
    list.dm_subtype <- list()
    for(subtype in subtypes){
        # SUBSET BY SUBTYPE ---
        dat2 <- dat1[which(dat1$subtype == subtype),]
        dat3 <- subset(dat2, select=c("motif","subtype","pval","sensitivity","fpr","ppv"))
        
        # PREPARE DATA FOR PLOT ---
            list.dm_subtype[[subtype]] <- prepareData(df=dat3, motif, subtype)

        cat("PROCESSED:", motif, subtype, "\n", sep=" ")        
    }

    # AGGREGATE DATA ---
    dm_subtype <- do.call(rbind.data.frame, list.dm_subtype)
    rownames(dm_subtype) <- NULL

    # FACTORIZE ---
    dm_subtype$pval <- factor(dm_subtype$pval, levels=c("p0",names(pvals),"p11"))
    dm_subtype$subtype <- factor(dm_subtype$subtype, levels=subtypes)

    # STORE DATA ---
    list.dm[[motif]] <- dm_subtype
}




### GENERATE PLOT ---
list.plot <- list()
for(motif in motifs){
    list.plot[[motif]] <- getROC(dm=list.dm[[motif]], motif, color_palette, text_label=FALSE) 
}



### GENERATE PLOT ---
file.plot <- file.path(dir.reproduce_fig, "supplementary_figure_13_a_c.pdf")
pdf(file.plot, height=2, width=6)
    grid.arrange(list.plot[[1]], list.plot[[2]], list.plot[[3]], nrow=1, ncol=3)    
dev.off()




##################################################################################################
#### CALCULATE AUC --------
list.res <- list()
for(motif in motifs){
    dm <- list.dm[[motif]]
    
    # GET AUC ---
    list.res[[motif]] <- getAUC(dm, subtypes, motif)
}

### GENERATE PLOT ---
list.plot_bar <- list()
for(motif in motifs){
    list.plot_bar[[motif]] <- plotBar(df=list.res[[motif]], subtypes, color_palette, motif)
    cat("PLOT GENERATED:", motif, "\n", sep=" ")     
}

### GENERATE PLOT ---
file.plot <- file.path(dir.reproduce_fig, "supplementary_figure_13_d_f.pdf")
pdf(file.plot, height=2, width=6)
    grid.arrange(list.plot_bar[[1]], list.plot_bar[[2]], list.plot_bar[[3]], nrow=1, ncol=3)
dev.off()


