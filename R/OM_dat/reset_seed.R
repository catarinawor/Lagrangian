#========================================================================================
# Generate input data for Lagrangian_gtg_OM
#========================================================================================





reset_seed<-function(DIR= "../../admb/OM/gtg/", sub_num=10){
	#"/Users/catarinawor/Documents/Lagrangian/admb/OM/gtg/"
	#setwd(DIR)
	myseed<-scan(file = "seed.txt")
	write.table(myseed-sub_num, file="seed.txt", row.names = F,
            col.names = F,)


}










