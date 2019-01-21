
library(shiny)
library(shinythemes)
library(DT)
putCap <- function(x, max, min){
  x[x > max] = max
  x[x < min] = min
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
                 
                 tags$h3("The following populations parameters are being simulated:"),
                 tags$h5(tags$b("Number of sections:"), "3"),
                 tags$h5(tags$b("Starting Population Size:"), "300"),
                 tags$h5(tags$b("Starting population fitness:"), "0.51"),
                 tags$h5(tags$b("Maximum population size during each day:"), "2000"),
                 tags$h5(tags$b("Survival threshold (environmental stress:"), "0.51"),
                 tags$h5(tags$b("Population mutation rate:"), "0.001"),
                 tags$br(),
                 tags$h3("Likelihood of survival:"),
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
                 
                 mainPanel(tabsetPanel(
                   tabPanel("Fitness Plot", plotOutput("sect1")),
                   tabPanel("Statistics and Summary", 
                            DT::dataTableOutput("tabFit.1"),
                            DT::dataTableOutput("predictFit.1"))
               )))
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
                 
                 mainPanel(tabsetPanel(
                   tabPanel("Fitness Plot", plotOutput("sect2")),
                   tabPanel("Statistics and Summary", 
                            DT::dataTableOutput("tabFit.2"),
                            DT::dataTableOutput("predictFit.2"))
                 )))
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
                 
                 mainPanel(tabsetPanel(
                   tabPanel("Fitness Plot", plotOutput("sect3")),
                   tabPanel("Statistics and Summary", 
                            DT::dataTableOutput("tabFit.3"),
                            DT::dataTableOutput("predictFit.3"))
                 )))
             )
    ) 
  ), 
  tabPanel(
    "Generations within Day",
    fluidPage(
      sidebarLayout(
        sidebarPanel(
          numericInput("Ng1.g", "Number of genes in first section",
                       value = 7),
          numericInput("Ng3.g", "Number of genes in third section",
                       value = 12),
          numericInput("Nl2.g", "Number of mutation sites in second-section genes",
                       value = 8),
          numericInput("Nl3.g", "Number of mutation sites in third-section genes",
                       value = 10),
          numericInput("day.g", "Days of Exposure to Antibiotics",
                      value = 5),
          numericInput("genpd", "Number of Generations per Day",
                        value = 24),
          sliderInput("ylim.g", label = "y-axis range",
                      min = 0,
                      max = 3000,
                      value = c(0, 2000),
                      step = 100)
          ),
          mainPanel(plotOutput("PopSize"))
      ))
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
  
  # Data table panel displays 
  output$tabFit.1 <- DT::renderDataTable({
    format.data.frame(data.frame(summary(d.lm1())$coefficients), digits = 4)
  })
  output$tabFit.2 <- DT::renderDataTable({
    format.data.frame(data.frame(summary(d.lm2())$coefficients), digits = 4)
  })
  output$tabFit.3 <- DT::renderDataTable({
    format.data.frame(data.frame(summary(d.lm3())$coefficients), digits = 4)
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
    putCap(predict(d.lm1(), newdata = newD()), 1,0)
  })
  fit.pred2 <- reactive({
    putCap(predict(d.lm2(), newdata = newD()), 1,0)
  })
  fit.pred3 <- reactive({
    putCap(predict(d.lm3(), newdata = newD()), 1,0)
  })
  fit.pred <- reactive({
    (fit.pred1() + fit.pred2() + fit.pred3())/3
  })
  
  
  # predict fitness based on user input for sect1 
  fit1.pred1 <- reactive({
    putCap(predict(d.lm1(), newdata = newD1.1()), 1,0)
  })
  
  output$predictFit.1 <- DT::renderDataTable({
    format.data.frame(data.frame(fit.pred1()), digits = 4)
  })
  
  # predict fitness based on user input for sect2 
  fit2.pred2 <- reactive({
    putCap(predict(d.lm2(), newdata = newD2.2()), 1,0)
  })
  
  output$predictFit.2 <- DT::renderDataTable({
    format.data.frame(data.frame(fit.pred2()), digits = 4)
  })
  
  # predict fitness based on user input for sect3 
  fit3.pred3 <- reactive({
    putCap(predict(d.lm3(), newdata = newD3.3()), 1,0)
  })
  
  output$predictFit.3 <- DT::renderDataTable({
    format.data.frame(data.frame(fit.pred3()), digits = 4)
  })
  
  
  success.pred <- reactive({
    predict(success.lm(), newdata = newD.success)
  })
  
  
  output$sect1 <- renderPlot({
    plot(1:input$days.1, fit1.pred1(), xlim = c(1, input$days.1), ylim = input$ylim.1, col='red', pch=16, cex = 2, cex.lab = 1.7, cex.axis = 1.5, 
         ylab = 'Section 1 Fitness', xlab = 'Days', type = "b", lty = 2)
  })
  
  output$sect2 <- renderPlot({
    plot(1:input$days.2, fit2.pred2(), xlim = c(1, input$days.2), ylim = input$ylim.2, col='dark orange', pch=16, cex = 2, cex.lab = 1.7, cex.axis = 1.5, 
         ylab = 'Section 2 Fitness', xlab = 'Days', type = "b", lty = 2)
  })
  
  output$sect3 <- renderPlot({
    plot(1:input$days.3, fit3.pred3(), xlim = c(1, input$days.3), ylim = input$ylim.3, col='yellow', pch=16, cex = 2, cex.lab = 1.7, cex.axis = 1.5, 
         ylab = 'Section 3 Fitness', xlab = 'Days', type = "b", lty = 2)
  })
  
  output$simulate <- renderPlot({
    par(oma = c(4, 1, 1, 1), mar = c(6, 5, 4, 2.5))
    
    plot(1:input$days, fit.pred1(), xlim = c(1, input$days), ylim = input$ylim, col='red', pch=16, cex = 2, cex.lab = 1.7, cex.axis = 1.5, 
         ylab = 'Overall and Section Fitness', xlab = 'Days', type = "b", lty = 2)
    points(1:input$days, fit.pred2(), col = "dark orange", pch=16, cex = 1.5, type = "b", lty = 2)
    points(1:input$days, fit.pred3(), col = "yellow", pch=16, cex = 1.5, type = "b", lty = 2)
    points(1:input$days, fit.pred(), col = "dark green", pch=16, cex = 1.5, type = "b", lty = 2)
    legend("bottom", inset = c(0, -0.65), legend=c("Overall Fitness", "Section 1 Fitness", "Section 2 Fitness", "Section 3 Fitness"),
           col=c("dark green", "red", "orange", "yellow"), pch=16, cex = 1.1, xpd = NA, ncol = 2, text.width = c(8,8,8,8), x.intersp = .2)
  })
  
  output$successRate <- renderText({
    paste("These population prameters are estimated to yield a population with a success rate of", predict(lm(success ~ Nl2 + Nl3 + Ng1 + Ng3, d), newdata = data.frame(Nl2 = input$Nl2, Nl3 = input$Nl3, Ng1 = input$Ng1, Ng3 = input$Ng3)), "%. Success rate is indicative of such a population's probability of survival.")
  })
  
## Generations in Each Day Plot
  fitGen.lm <- reactive({
    lm(fit ~ poly(gen.number,3, raw = T) + Nl2 + Nl3 + Ng1 + Ng3 + Day, genDay[genDay$Day == input$day.g,])
  })    
  dataGenFit <- reactive({
    data.frame(Day = input$day.g, Nl2 = input$Nl2, Nl3 = input$Nl3, Ng1 = input$Ng1, Ng3 = input$Ng3, gen.number = 1:input$genpd)
  })
  fitGen.pred <- reactive({
    predict(fitGen.lm(), newdata = dataGenFit())
  })
  
  genDay.lm <- reactive({
    lm(ni ~ poly(gen.number,3, raw = T) + Nl2 + Nl3 + Ng1 + Ng3 + fit + Day, genDay[genDay$Day == input$day.g,])
  })
  dataGenDay <- reactive({
    data.frame(Day = input$day.g, Nl2 = input$Nl2, Nl3 = input$Nl3, Ng1 = input$Ng1, Ng3 = input$Ng3, gen.number = 1:input$genpd, fit = fitGen.pred())
  })
  genDay.pred <- reactive({
    predict(genDay.lm(), newdata = dataGenDay())
  })
  
  output$PopSize <- renderPlot({
    plot(1:input$genpd, putCap(genDay.pred(), 2000, 0), ylim = input$ylim.g, col='navy blue', pch=16, cex = 2, cex.lab = 1.7, cex.axis = 1.5, 
         ylab = 'Population Size', xlab = 'Generation Number', type = "b", lty = 2)
  }) 
}

# Run the application
shinyApp(ui = ui, server = server)
