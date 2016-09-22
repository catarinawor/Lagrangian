#============================================================
# Bunch of useful functions stolen from various places
#============================================================


#Stolen for Cgradin
install.packages.if.needed <- function(package.name, package.install.name, github=FALSE){
  if(github){
    if(!(package.name %in% rownames(installed.packages()))){
      devtools::install_github(package.install.name)
    }
  }else{
    if(!(package.name %in% rownames(installed.packages()))){
      install.packages(package.install.name)
    }
  }
}