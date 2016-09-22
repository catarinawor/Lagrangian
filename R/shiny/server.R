source('helpers.R')


shinyServer(function(input,output,session){

    ## ------------------------------------------------------------ ##
    ## PACIFIC HAKE MIGRATION MODEL (NOV 10, 2015)
    ## ------------------------------------------------------------ ##
  
  #getLagrA<-reactive(do.call(getResult,getArgsLagr(input,"A")))
  #getLagrB<-reactive(do.call(getResult,getArgsLagr(input,"B")))
  
  getLagrA<-eventReactive(input$runButton,do.call(getResult,getArgsLagr(input,"A")))
  getLagrB<-eventReactive(input$runButton,do.call(getResult,getArgsLagr(input,"B")))
  


  
  print("cheguei aqui")

  output$lagr_plot<-renderPlot({ 
    
        MresA<-getLagrA()
        MresB<-getLagrB()
        plotvb(MresA,MresB)

    })

  output$catch_plot<-renderPlot({ 

        MresA<-getLagrA()
        MresB<-getLagrB()
        plotcatch(MresA,MresB)

    })
  





})
# End of shinyServer
