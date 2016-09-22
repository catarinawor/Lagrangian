## RENDER LAGRAGIAN USER INTERFACE
getArgsLagr <- function(input,prefix){

    print("in getargs")

    argsLagr <- list(sl_effm=input[[paste0(prefix,"_","sl_effm")]],rb_envr=input[[paste0(prefix,"_","rb_envr")]])
    print(argsLagr)
    
    return(argsLagr)
}


getResult<-function(sl_effm,rb_envr){

    
    print("in getResult")

    MP$effm<<-sl_effm

    if(rb_envr=="warm"){
        MP$maxPos50 <<- 3
        
    }else if(rb_envr=="cold"){
        
        MP$maxPos50 <<- 6
    
    }else{
        
        MP$maxPos50 <<- 4
    } 

    out<-lagr_func(MP)     
    
    return(out)
}