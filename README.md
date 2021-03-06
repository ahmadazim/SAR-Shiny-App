
# Shiny App Visualization for a Stochastic Simulation of Antibiotic Resistance [S.A.R.]

This project is a deterministic version of the **SAR** simulation package by creating a user-friendly interface for users to immediately visualize the effects of varying population parameters on the fitness of a population when placed under environmental stress. The SAR simulation package was originally written to answer questions addressed in the manuscript titled *Genetic basis and patterns of gaining antibiotic resistance of a naïve strain of Staphylococcus aureus*, authored by S., A., and G. Abdel-Azim, and is available at https://github.com/SalmaGAA/SAR. The study analyzed the pattern of amoxicillin-resistance in *Staphylococcus aureus* and, through the use of simulation, offered explanations for the observed step-wise pattern. The simulation package was useful in providing evidence of this step-wise pattern caused by varying levels of complexity in the mutations that bacteria undergo to gain resistance. The ability to set unique initial bacterial population parameters provides flexibility in its usage; however, it lacks a simple user-friendly interface, providing impetus for this Shiny app. The Shiny app is available on shinyapp.io at https://ahmadazim.shinyapps.io/shinySAR/.

# To download and launch Shiny app on your local system:
**Method 1:** download and run using *runGitHub()*
```r
library(shiny)
runGitHub("ahmadazim/SAR-Shiny-App", "rstudio")
```

**Method 2:** download the source code and run the app\
First download the git repository (for example, to "~/SAR-Shiny-App").
```r
# Run the following in RStudio console:
setwd("~/SAR-Shiny-App")
shiny::runApp()
```
