## RENDER LAGRAGIAN USER INTERFACE
renderLAGR <- function(){	
	   buildLAGRGui()
       
}



buildLAGRGui  <- function(){

   fluidPage( 
       
        column(5,

            #buildTMAInterface()   
            wellPanel( 
                buildLagrSCN(),
                actionButton("runButton", "Run!")    
            )
        ),
            
        column(7,                                                
            wellPanel( 
                plotOutput("lagr_plot"), 
                plotOutput("catch_plot")
            )

            
               
                         
        )                   
   ) 
                
}

buildLagrInputs<-function(prefix){
    fluidRow(
                    

        wellPanel("Fishing Effort Intensity",                  

            sliderInput(paste0(prefix,"_","sl_effm"), "Effort Intensity", 1,10,1)        
        ),

        wellPanel("Environmental regime", 
            radioButtons(paste0(prefix,"_","rb_envr"), "Choose environmental option:",
            c("warm", "average","cold"),selected = "average")
        )


    )     
}

buildLagrSCN <- function(){
    
    fluidRow(
    
        column(6,
            tags$h4("Scenario A"),
            buildLagrInputs("A")
            ), 

        column(6,
            tags$h4("Scenario B"),
            buildLagrInputs("B")
            )       
    )

}


