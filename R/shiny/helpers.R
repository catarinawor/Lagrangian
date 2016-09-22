# 
# Load Libraries
# 
if(!require("shiny"))         install.packages("shiny")
if(!require("shinythemes"))   install.packages("shinythemes")
if(!require("ggplot2"))       install.packages("ggplot2")
if(!require("reshape2"))      install.packages("reshape2")
if(!require("Rcpp"))          install.packages("Rcpp")
if(!require("markdown"))      install.packages("markdown")
if(!require("plyr"))          install.packages("plyr")
if(!require("dplyr"))         install.packages("dplyr")
if(!require("grid"))          install.packages("grid")
if(!require("shinydashboard"))install.packages("shinydashboard")
if(!require("dygraphs"))      install.packages("dygraphs")
if(!require("devtools"))          install.packages("devtools")
if(!require("shinyTable"))    install_github("shinyTable", "trestletech")



# Source R-scripts.

.LIB        <- "./data"
.RFILES     <- list.files(.LIB,pattern="\\.[Rr]$",recursive=TRUE,full.names=TRUE)
for(nm in .RFILES){
        print(nm)
        source(file.path(nm), echo=FALSE, local=TRUE)
}        


#
# Source Rcpp scripts
# 

