cat("###############################################################################################################\n", sep=" ")
cat(format(Sys.time(), "%b %d %X"), "LOADAING PATHS ... ", "\n", sep=" ")

#### DEFINE PATHS: PATH WHERE ATAC-SEQ BAMs/PEAK FILES ARE STORED ------------------------
dir.base <- file.path("/data1/datasets_1/human_prostate_WCDT/atac_solid")
dir.process <- file.path(dir.base, "processed/2020_09_22/processed") # CHANGE THIS AS REQUIRED: PATH TO THE PROCESSED SAMPLE FOLDER CONTAINING BAM AND PEAK FILES

#### DEFINE PATHS: PATH DATA ANALYSIS PROJECT --------------------------------------------
dir.wrk <- file.path("/data1/projects/WCDT_atac_2020")
dir.metadata <- file.path(dir.wrk, "metadata")
dir.reproduce <- file.path(dir.wrk, "reproduce") # CHANGE THIS AS REQUIRED: PATH TO RE-PRODUCTION DIRECTORY
dir.reproduce_data <- file.path(dir.reproduce, "data") # PATH TO DIRECTORY CONTAINING ALL REQUIRED DATA 
dir.reproduce_fig <- file.path(dir.reproduce, "figures") # PATH TO DIRECTORY CONTAINING ALL OUTPUT MAIN AND SUPPLEMENTARY FIGURES
dir.reproduce_tbl <- file.path(dir.reproduce, "tables") # PATH TO DIRECTORY CONTAINING ALL OUTPUT SUPPLEMENTARY TABLES
dir.reproduce_scripts <- file.path(dir.reproduce, "scripts")  # PATH TO DIRECTORY CONTAINING SCRIPTS

#### SET WORKING DIRECTORY TO PRODUCTION DIRECTORY ---------------------------------------
setwd(dir.reproduce)
cat(format(Sys.time(), "%b %d %X"), "WORKING DIRECTORY IS SET TO:", getwd(), "\n", sep=" ")
cat(format(Sys.time(), "%b %d %X"), "dir.reproduce_data:", dir.reproduce_data, "\n", sep=" ")
cat(format(Sys.time(), "%b %d %X"), "dir.reproduce_fig:", dir.reproduce_fig, "\n", sep=" ")
cat(format(Sys.time(), "%b %d %X"), "dir.reproduce_tbl:", dir.reproduce_tbl, "\n", sep=" ")
cat(format(Sys.time(), "%b %d %X"), "dir.reproduce_scripts:", dir.reproduce_scripts, "\n", sep=" ")

cat(format(Sys.time(), "%b %d %X"), "DONE!", "\n", sep=" ")
cat("###############################################################################################################\n", sep=" ")
