
library(shiny)
library(DT)
library(SAR)
library(shinythemes)
putCap <- function(x){
  x[x > 1] = 1
  x[x < 0] = 0
  return(x)
}


ui <- navbarPage(
  theme = shinytheme("sandstone"),
  "Effect of Bacterial Genotype on Growth Levels Over Multiple Exposures",
  tabPanel("Simulation",
  fluidPage(
   sidebarLayout(
      sidebarPanel(
        numericInput("Ng1", "Number of genes in first section",
                     value = 1),
        numericInput("Ng2", "Number of genes in second section",
                     value = 7),
        numericInput("Ng3", "Number of genes in third section",
                     value = 15),
        numericInput("Nl1", "Number of mutation sites in first-section genes",
                     value = 1),
        numericInput("Nl2", "Number of mutation sites in second-section genes",
                     value = 6),
        numericInput("Nl3", "Number of mutation sites in third-section genes",
                     value = 11),
        sliderInput("days", "Length of Simulation, in days",
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
        textOutput("parameters"),
        textOutput("sect.text"),
        textOutput("Nl1.text"),
        textOutput("Nl2.text"),
        textOutput("Nl3.text"),
        textOutput("Ng1.text"),
        textOutput("Ng2.text"),
        textOutput("Ng3.text"),
        textOutput("startPop"),
        textOutput("startPopFit"),
        textOutput("maxPop"),
        textOutput("nDays"),
        textOutput("thr"),
        textOutput("mr"),
        textOutput("successRate")
        )
   )
  )
  ),
  navbarMenu(
    "Fitness Levels by Section",
      tabPanel("Section 1 Fitness Levels",
             fluidPage(
               sidebarLayout(
                 sidebarPanel(
                   numericInput("Ng1.1", "Number of genes in first section",
                                value = 7),
                   numericInput("Ng2.1", "Number of genes in second section",
                                value = 10),
                   numericInput("Ng3.1", "Number of genes in third section",
                                value = 12),
                   numericInput("Nl1.1", "Number of mutation sites in first-section genes",
                                value = 6),
                   numericInput("Nl2.1", "Number of mutation sites in second-section genes",
                                value = 8),
                   numericInput("Nl3.1", "Number of mutation sites in third-section genes",
                                value = 10),
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
                   numericInput("Ng1.2", "Number of genes in first section",
                                value = 7),
                   numericInput("Ng2.2", "Number of genes in second section",
                                value = 10),
                   numericInput("Ng3.2", "Number of genes in third section",
                                value = 12),
                   numericInput("Nl1.2", "Number of mutation sites in first-section genes",
                                value = 6),
                   numericInput("Nl2.2", "Number of mutation sites in second-section genes",
                                value = 8),
                   numericInput("Nl3.2", "Number of mutation sites in third-section genes",
                                value = 10),
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
                   numericInput("Ng1.3", "Number of genes in first section",
                                value = 7),
                   numericInput("Ng2.3", "Number of genes in second section",
                                value = 10),
                   numericInput("Ng3.3", "Number of genes in third section",
                                value = 12),
                   numericInput("Nl1.3", "Number of mutation sites in first-section genes",
                                value = 6),
                   numericInput("Nl2.3", "Number of mutation sites in second-section genes",
                                value = 8),
                   numericInput("Nl3.3", "Number of mutation sites in third-section genes",
                                value = 10),
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
     
  d.lm1 <- reactive({
    lm(fit.1 ~ poly(days,3, raw = T) + Nl1 + Nl2 + Nl3 + Ng1 + Ng2 + Ng3, d)
    })
  d.lm2 <- reactive({
    lm(fit.2 ~ poly(days,3, raw = T) + Nl1 + Nl2 + Nl3 + Ng1 + Ng2 + Ng3, d)
    })
  d.lm3 <- reactive({
    lm(fit.3 ~ poly(days,3, raw = T) + Nl1 + Nl2 + Nl3 + Ng1 + Ng2 + Ng3, d)
    })
  
  success.lm <- reactive({
    lm(success ~ + Nl1 + Nl2 + Nl3 + Ng1 + Ng2 + Ng3, d)
     })
  
  # Construct new data using user input for overall
  newD <- reactive({
    data.frame(days = 1:input$days, Nl1 = input$Nl1, Nl2 = input$Nl2, Nl3 = input$Nl3, Ng1 = input$Ng1, Ng2 = input$Ng2, Ng3 = input$Ng3)
  })
  
  # construct new data using user input for sections
  newD1.1 <- reactive({
    data.frame(days = 1:input$days.1, Nl1 = input$Nl1.1, Nl2 = input$Nl2.1, Nl3 = input$Nl3.1, Ng1 = input$Ng1.1, Ng2 = input$Ng2.1, Ng3 = input$Ng3.1)
  })
  newD2.2 <- reactive({
    data.frame(days = 1:input$days.2, Nl1 = input$Nl1.2, Nl2 = input$Nl2.2, Nl3 = input$Nl3.2, Ng2 = input$Ng2.2, Ng1 = input$Ng1.2, Ng3 = input$Ng3.2)
  })
  newD3.3 <- reactive({
    data.frame(days = 1:input$days.3, Nl1 = input$Nl1.3, Nl2 = input$Nl2.3, Nl3 = input$Nl3.3, Ng1 = input$Ng1.3, Ng2 = input$Ng2.3, Ng3 = input$Ng3.3)
  })
  
  newD.success <- reactive({
    data.frame(Nl1 = input$Nl1, Nl2 = input$Nl2, Nl3 = input$Nl3, Ng1 = input$Ng1, Ng2 = input$Ng2, Ng3 = input$Ng3)
  })
  
  
  
  # predict fitness based on user input for overall 
  fit.pred1 <- reactive({
    putCap(predict(d.lm1(), newdata = newD()))
  })
  fit.pred2 <- reactive({
    putCap(predict(d.lm2(), newdata = newD()))
  })
  fit.pred3 <- reactive({
    putCap(predict(d.lm3(), newdata = newD()))
  })
  fit.pred <- reactive({
    (fit.pred1() + fit.pred2() + fit.pred3())/3
  })
  
  
  # predict fitness based on user input for sect1 
  fit1.pred1 <- reactive({
    putCap(predict(d.lm1(), newdata = newD1.1()))
  })
  
  # predict fitness based on user input for sect2 
  fit2.pred2 <- reactive({
    putCap(predict(d.lm2(), newdata = newD2.2()))
  })
  
  # predict fitness based on user input for sect3 
  fit3.pred3 <- reactive({
    putCap(predict(d.lm3(), newdata = newD3.3()))
  })
  
  
  
  success.pred <- reactive({
    predict(success.lm(), newdata = newD.success)
  })
  
  
  output$sect1 <- renderPlot({
    plot(1:input$days.1, fit1.pred1(), xlim = c(1, input$days.1), ylim = input$ylim.1, col='red', pch=16, cex = 2, cex.lab = 1.7, cex.axis = 1.5, 
         ylab = 'Section 1 Fitness', xlab = 'Days', type = "b")
  })
  
  output$sect2 <- renderPlot({
    plot(1:input$days.2, fit2.pred2(), xlim = c(1, input$days.2), ylim = input$ylim.2, col='dark orange', pch=16, cex = 2, cex.lab = 1.7, cex.axis = 1.5, 
         ylab = 'Section 2 Fitness', xlab = 'Days', type = "b")
  })
  
  output$sect3 <- renderPlot({
    plot(1:input$days.3, fit3.pred3(), xlim = c(1, input$days.3), ylim = input$ylim.3, col='yellow', pch=16, cex = 2, cex.lab = 1.7, cex.axis = 1.5, 
         ylab = 'Section 3 Fitness', xlab = 'Days', type = "b")
  })
  
  output$simulate <- renderPlot({
    par(oma = c(4, 1, 1, 1))
    
    plot(1:input$days, fit.pred1(), xlim = c(1, input$days), ylim = input$ylim, col='red', pch=16, cex = 2, cex.lab = 1.7, cex.axis = 1.5, 
         ylab = 'Overall and Section Fitness', xlab = 'Days', type = "b")
    points(1:input$days, fit.pred2(), col = "dark orange", pch=16, cex = 1.5)
    points(1:input$days, fit.pred3(), col = "yellow", pch=16, cex = 1.5)
    points(1:input$days, fit.pred(), col = "dark green", pch=16, cex = 1.5)
    legend("bottom", inset = c(0, -0.55), legend=c("Overall Fitness", "Section 1 Fitness", "Section 2 Fitness", "Section 3 Fitness"),
           col=c("dark green", "red", "orange", "yellow"), pch=16, cex = 1.1, xpd = NA, ncol = 2, text.width = c(8,8,8,8), x.intersp = .2)
  })
  
  output$parameters <- renderText({
    paste("The following populations parameters are being simulated...")
  })

  output$sect.text <- renderText({
    paste("Number of sections: 3")
  })
  
  output$Nl1.text <- renderText({
    paste("Number of mutation sites in the first-section genes:", input$Nl1)
  })
  
  output$Nl2.text <- renderText({
    paste("Number of mutation sites in the second-section genes:", input$Nl2)
  })
  
  output$Nl3.text <- renderText({
    paste("Number of mutation sites in the third-section genes:", input$Nl3)
  })
  
  output$Ng1.text <- renderText({
    paste("Number of genes in the first section:", input$Ng1)
  })
  
  output$Ng2.text <- renderText({
    paste("Number of genes in the second section:", input$Ng2)
  })
  
  output$Ng3.text <- renderText({
    paste("Number of genes in the third section:", input$Ng3)
  })
  
  output$startPop <- renderText({
    paste("Starting Population Size: 300")
  })
  
  output$startPopFit <- renderText({
    paste("Starting population fitness: 0.51")
  })
  
  output$maxPop <- renderText({
    paste("Maximum population size during each day: 2000")
  })
  
  output$nDays <- renderText({
    paste("Number of days being simulated:", input$days)
  })
  
  output$thr <- renderText({
    paste("Survival threshold (environmental stress): 0.51")
  })

  output$mr <- renderText({
    paste("Population mutation rate: 0.001")
  })
  
  output$successRate <- renderText({
    paste("These population prameters are estimated to yield a population with a success rate of", predict(lm(success ~ Nl2 + Nl3 + Ng1 + Ng3, d), newdata = data.frame(Nl2 = input$Nl2, Nl3 = input$Nl3, Ng1 = input$Ng1, Ng3 = input$Ng3)), "%. Success rate is indicative of such a population's probability of survival.")
  })
  
}

# Run the application
shinyApp(ui = ui, server = server)
