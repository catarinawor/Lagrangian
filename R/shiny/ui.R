# ----------------------------------------#
# MAIN USER INTERFACE FOR THE APPLICATION #
# ----------------------------------------#
source("helpers.R")
shinyUI(fluidPage(
        navbarPage(
            "Pacific hake migration model",

            ## -------------- ##
            # About Interface  #
            ## -------------- ##
            #renderAbout(),

            ## -------------- ##
            # Equil Interface  #
            ## -------------- ##
            renderLAGR()

            ## -------------- ##
            # MAPS Interface   #
            ## -------------- ##
            #renderMaps(),

        )
    )
)
