#==============================================
#Title: Read in simulation outcomes and plot thes
#Author: Catarina Wor
#date: Oct 6 2015
#
#==============================================


library(ggplot2)
library(reshape2)



.PWD        <- "/Users/catarinawor/Documents/Lagrangian"
.THEME      <- theme_bw(11)

.SIMDIRS   <- "/Users/catarinawor/Documents/Lagrangian/SimResult"
.SIMNAME   <- list.files(.SIMDIRS,pattern="\\.Rdata", full.name=TRUE)
load(.SIMNAME)

load("SimResult/simest5262.Rdata")
names(M)     <- strsplit(basename(.SIMNAME),".RData")



?list.files


