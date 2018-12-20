
library(shiny)
library(DT)
library(SAR)
library(shinythemes)


ui <- navbarPage(
  theme = shinytheme("sandstone"),
  "Effect of Bacterial Genotype on Growth Levels Over Multiple Exposures",
  tabPanel("Simulation",
  fluidPage(
   sidebarLayout(
      sidebarPanel(
        numericInput("Nl2", "Number of loci in second section",
                     value = 6),
        numericInput("Nl3", "Number of loci in third section",
                     value = 10),
        numericInput("Ng1", "Number of genes in first-section loci",
                    value = 7),
        numericInput("Ng3", "Number of genes in third-section loci",
                    value = 12),
        sliderInput("days", "Length of Simulation, in Days",
                    min = 2,
                    max = 40,
                    value = 30),
        sliderInput("ylim", label = "y-axis range",
                    min = 0,
                    max = 3,
                    value = c(0, 1.2),
                    step = 0.1)
      ),

      mainPanel(
        plotOutput("simulate"),
        textOutput("successRate")
        )
   )
  )
  ),
  navbarMenu(
    "Fitness Levels by Section",
    tabPanel("Overall Genome Fitness Levels",
             fluidPage(
               sidebarLayout(
                 sidebarPanel(
                   numericInput("Nl2.O", "Number of loci in second section",
                                value = 6),
                   numericInput("Nl3.O", "Number of loci in third section",
                                value = 10),
                   numericInput("Ng1.O", "Number of genes in first-section loci",
                                value = 7),
                   numericInput("Ng3.O", "Number of genes in third-section loci",
                                value = 12),
                   sliderInput("days.O", "Length of Simulation, in Days",
                               min = 2,
                               max = 40,
                               value = 30),
                   sliderInput("ylim.O", label = "y-axis range",
                               min = 0,
                               max = 3,
                               value = c(0, 1.2),
                               step = 0.1)
                 ),
                 
                 mainPanel(plotOutput("overall"))
               )
              )
    ),
    tabPanel("Section 1 Fitness Levels",
             fluidPage(
               sidebarLayout(
                 sidebarPanel(
                   numericInput("Nl2.1", "Number of loci in second section",
                                value = 6),
                   numericInput("Nl3.1", "Number of loci in third section",
                                value = 10),
                   numericInput("Ng1.1", "Number of genes in first-section loci",
                                value = 7),
                   numericInput("Ng3.1", "Number of genes in third-section loci",
                                value = 12),
                   sliderInput("days.1", "Length of Simulation, in Days",
                               min = 2,
                               max = 40,
                               value = 30),
                   sliderInput("ylim.1", label = "y-axis range",
                               min = 0,
                               max = 3,
                               value = c(0, 1.2),
                               step = 0.1)
                 ),
                 
                 mainPanel(plotOutput("sect1"))
               )
             )
    ),
    tabPanel("Section 2 Fitness Levels",
             fluidPage(
               sidebarLayout(
                 sidebarPanel(
                   numericInput("Nl2.2", "Number of loci in second section",
                                value = 6),
                   numericInput("Nl3.2", "Number of loci in third section",
                                value = 10),
                   numericInput("Ng1.2", "Number of genes in first-section loci",
                                value = 7),
                   numericInput("Ng3.2", "Number of genes in third-section loci",
                                value = 12),
                   sliderInput("days.2", "Length of Simulation, in Days",
                               min = 2,
                               max = 40,
                               value = 30),
                   sliderInput("ylim.2", label = "y-axis range",
                               min = 0,
                               max = 3,
                               value = c(0, 1.2),
                               step = 0.1)
                 ),
                 
                 mainPanel(plotOutput("sect2"))
               )
             )
    ),
    tabPanel("Section 3 Fitness Levels",
             fluidPage(
               sidebarLayout(
                 sidebarPanel(
                   numericInput("Nl2.3", "Number of loci in second section",
                                value = 6),
                   numericInput("Nl3.3", "Number of loci in third section",
                                value = 10),
                   numericInput("Ng1.3", "Number of genes in first-section loci",
                                value = 7),
                   numericInput("Ng3.3", "Number of genes in third-section loci",
                                value = 12),
                   sliderInput("days.3", "Length of Simulation, in Days",
                               min = 2,
                               max = 40,
                               value = 30),
                   sliderInput("ylim.3", label = "y-axis range",
                               min = 0,
                               max = 3,
                               value = c(0, 1.2),
                               step = 0.1)
                 ),
                 
                  mainPanel(plotOutput("sect3"))
               )
             )
    ) 
  ), collapsible = TRUE
)
  

server <- function(input, output) {
     
  d.lm <- reactive({
    lm(fit ~ poly(days,3, raw = T) + Nl2 + Nl3 + Ng1 + Ng3, d)
  })
  d.lm1 <- reactive({
    lm(fit.1 ~ poly(days,3, raw = T) + Nl2 + Nl3 + Ng1 + Ng3, d)
    })
  d.lm2 <- reactive({
    lm(fit.2 ~ poly(days,3, raw = T) + Nl2 + Nl3 + Ng1 + Ng3, d)
    })
  d.lm3 <- reactive({
    lm(fit.3 ~ poly(days,3, raw = T) + Nl2 + Nl3 + Ng1 + Ng3, d)
    })
  
  success.lm <- reactive({
    lm(success ~ Nl2 + Nl3 + Ng1 + Ng3, d)
     })
  
  # construct new data using user input
  newD <- reactive({
    data.frame(days = 1:input$days, Nl2 = input$Nl2, Nl3 = input$Nl3, Ng1 = input$Ng1, Ng3 = input$Ng3)
    })
  newDO <- reactive({
    data.frame(days = 1:input$days.O, Nl2 = input$Nl2.O, Nl3 = input$Nl3.O, Ng1 = input$Ng1.O, Ng3 = input$Ng3.O)
  })
  newD1 <- reactive({
    data.frame(days = 1:input$days.1, Nl2 = input$Nl2.1, Nl3 = input$Nl3.1, Ng1 = input$Ng1.1, Ng3 = input$Ng3.1)
  })
  newD2 <- reactive({
    data.frame(days = 1:input$days.2, Nl2 = input$Nl2.2, Nl3 = input$Nl3.2, Ng1 = input$Ng1.2, Ng3 = input$Ng3.2)
  })
  newD3 <- reactive({
    data.frame(days = 1:input$days.3, Nl2 = input$Nl2.3, Nl3 = input$Nl3.3, Ng1 = input$Ng1.3, Ng3 = input$Ng3.3)
  })
  
  newD.success <- reactive({
    data.frame(Nl2 = input$Nl2, Nl3 = input$Nl3, Ng1 = input$Ng1, Ng3 = input$Ng3)
  })
  
  
  
  # predict fitness based on user input 
  fit.pred <- reactive({
    predict(d.lm(), newdata = newD())
    })
  fit.predO <- reactive({
    predict(d.lm(), newdata = newDO())
  })
  fit.pred1 <- reactive({
    predict(d.lm1(), newdata = newD1())
  })
  fit.pred2 <- reactive({
    predict(d.lm2(), newdata = newD2())
  })
  fit.pred3 <- reactive({
    predict(d.lm3(), newdata = newD3())
  })
  
  success.pred <- reactive({
    predict(success.lm(), newdata = newD.success)
  })
  
  
  output$overall <- renderPlot({
    plot(1:input$days.O, fit.predO(), xlim = c(1, input$days.O), ylim = input$ylim.O, col='dark green', pch=16, cex = 2, cex.lab = 1.7, cex.axis = 1.5, 
         ylab = 'Overall Fitness', xlab = 'Days', type = "b")
  })
  
  output$sect1 <- renderPlot({
    plot(1:input$days.1, fit.pred1(), xlim = c(1, input$days.1), ylim = input$ylim.1, col='red', pch=16, cex = 2, cex.lab = 1.7, cex.axis = 1.5, 
         ylab = 'Section 1 Fitness', xlab = 'Days', type = "b")
  })
  
  output$sect2 <- renderPlot({
    plot(1:input$days.2, fit.pred2(), xlim = c(1, input$days.2), ylim = input$ylim.2, col='dark orange', pch=16, cex = 2, cex.lab = 1.7, cex.axis = 1.5, 
         ylab = 'Section 2 Fitness', xlab = 'Days', type = "b")
  })
  
  output$sect3 <- renderPlot({
    plot(1:input$days.3, fit.pred3(), xlim = c(1, input$days.3), ylim = input$ylim.3, col='yellow', pch=16, cex = 2, cex.lab = 1.7, cex.axis = 1.5, 
         ylab = 'Section 3 Fitness', xlab = 'Days', type = "b")
  })
  
  output$simulate <- renderPlot({
    par(oma = c(4, 1, 1, 1))
    
    plot(1:input$days, fit.pred(), xlim = c(1, input$days), ylim = input$ylim, col='dark green', pch=16, cex = 2, cex.lab = 1.7, cex.axis = 1.5, 
         ylab = 'Overall and Section Fitness', xlab = 'Days', type = "b")
    points(1:input$days.1, fit.pred1(), col = "red", pch=16, cex = 1.5)
    points(1:input$days.2, fit.pred2(), col = "dark orange", pch=16, cex = 1.5)
    points(1:input$days.3, fit.pred3(), col = "yellow", pch=16, cex = 1.5)
    legend("bottom", inset = c(0, -0.55), legend=c("Overall Fitness", "Section 1 Fitness", "Section 2 Fitness", "Section 3 Fitness"),
           col=c("dark green", "red", "orange", "yellow"), pch=16, cex = 1.1, xpd = NA, ncol = 2, text.width = c(8,8,8,8), x.intersp = .2)
  })
  
  output$successRate <- renderText({
    paste("These population prameters are estimated to yield a population with a success rate of", predict(lm(success ~ Nl2 + Nl3 + Ng1 + Ng3, d), newdata = data.frame(Nl2 = input$Nl2, Nl3 = input$Nl3, Ng1 = input$Ng1, Ng3 = input$Ng3)), "%. Success rate is indicative of such a population's probability of survival.")
    })
}

# Run the application
shinyApp(ui = ui, server = server)
