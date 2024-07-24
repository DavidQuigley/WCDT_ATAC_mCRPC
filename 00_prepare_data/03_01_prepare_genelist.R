###################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2023 ##
## R. SHRESTHA, ET AL., D.A. QUIGLEY AND F.Y. FENG LABS, UCSF    ##
###################################################################

### DEFINE LIBRARIES ---
library("stringr")
library("data.table")
library("openxlsx")

### DEFINE PATH ---
dir.reproduce <- file.path("/data1/projects/WCDT_atac_2020/reproduce")
file.config_path <- file.path(dir.reproduce, "scripts/config/config_paths.R") 
source(file.config_path)

### DEFINE FILES ---
file.msigdb_biocarta <- file.path(dir.wrk, "reference/msigdb/msigdb_v7.5.1/c2.cp.biocarta.v7.5.1.symbols.gmt")
file.msigdb_kegg <- file.path(dir.wrk, "reference/msigdb/msigdb_v7.5.1/c2.cp.kegg.v7.5.1.symbols.gmt")
file.msigdb_reactome <- file.path(dir.wrk, "reference/msigdb/msigdb_v7.5.1/c2.cp.reactome.v7.5.1.symbols.gmt")
file.msigdb_wikipathways <- file.path(dir.wrk, "reference/msigdb/msigdb_v7.5.1/c2.cp.wikipathways.v7.5.1.symbols.gmt")
file.msigdb_gobp <- file.path(dir.wrk, "reference/msigdb/msigdb_v7.5.1/c5.go.bp.v7.5.1.symbols.gmt")
file.msigdb_hallmark <- file.path(dir.wrk, "reference/msigdb/msigdb_v7.5.1/h.all.v7.5.1.symbols.gmt")
file.nepc_data <- file.path(dir.wrk, "external_data/beltran_nepc_2016/supplementary/beltran_etal_nature_medicine_2016_supplementary_tables.xlsx")
file.tang2022 <- file.path(dir.wrk, "reference/genelist/genelist_tang_etal_science_markers.tsv")
file.labrecque2019 <- file.path(dir.wrk, "reference/genelist/genelist_nepc.RData")
file.utility_functions <- file.path(dir.reproduce_scripts, "00_prepare_data/00_01_utility_functions.R")

### LOAD FUNCTIONS ---
source(file.utility_functions)



##################################################################################################################
#### GET MSIGDB PATHWAYS ----------
msigdb.pathways <- c("biocarta","kegg","reactome","wikipathways","gobp","hallmark")
file.msigdb <- c(file.msigdb_biocarta, file.msigdb_kegg, file.msigdb_reactome, file.msigdb_wikipathways, file.msigdb_gobp, file.msigdb_hallmark)
names(file.msigdb) <- msigdb.pathways

list.msigdb_pathways <- list()
for(item in msigdb.pathways){
    list.msigdb_pathways[[item]] <- gmtPathways(file.gmt=file.msigdb[[item]])
}





##################################################################################################################
### GET NEPC GENES ---
genes_nepc <- openxlsx::read.xlsx(xlsxFile=file.nepc_data, sheet=10, cols=1, startRow=2)$HGNC.ID


### GET TANG et al 2022 Science GENELIST ---
dat_tang2022 <- data.table::fread(file=file.tang2022, sep="\t", header=TRUE, nThread=1,  data.table=FALSE, verbose=FALSE)
#dat_tang2022 <- subset(dat_tang2022, dat_tang2022$Category %in% c("AR","WNT","NE","Stem_cell"))

### GET Labrecque et al 2019 GENELIST ---
load(file.labrecque2019)
dat_labrecque2019 <- NEPC.gene.list$Nelson
colnames(dat_labrecque2019) <- c("Gene","Phenotype","Color")


###################################################################################################################
### ADD TO LIST ---
list.output <- list(msigdb_pathways=list.msigdb_pathways,
                    genes_nepc=genes_nepc,
                    genes_tang2022=dat_tang2022,
                    genes_labrecque2019=dat_labrecque2019)


### SAVE OBJECT TO RDATA FILE ---
file.rds_genelist <- file.path(dir.reproduce_data, "genelist.rds")
saveRDS(object=list.output, file=file.rds_genelist)
