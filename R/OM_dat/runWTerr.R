#========================================================================================
# Generate input data for Lagrangian_gtg_OM
#========================================================================================



args <- commandArgs(trailingOnly=TRUE)


#source("/Users/catarinawor/Documents/Lagrangian/R/OM_dat/rec_dat.R")
source("rec_dat.R")


sed<-scan("../../admb/OM/gtg/seed.txt")


writeWTdat(myseed=sed)











