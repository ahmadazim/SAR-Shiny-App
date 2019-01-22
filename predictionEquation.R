cleanUpZeros <- function(ar) {
  # ar: 3-dim'l array
  i0 <- which(apply(ar, 3, function(x) {mean(x[,1])}) == 0)
  if(length(i0)) return(ar[,,-i0])
  else return(ar)
}

## Flatten into 2-dim'l and characterize resulting array
flaten2D <- function(ar, Nl, Ng) {
  successRuns <- dim(ar)[3]
  TwoDimAr <- do.call('rbind', lapply(1:successRuns, function(x) ar[,,x]))
  ## attach number of successful runs, number of loci and number of genes to the 2-dim'l flat array
  TwoDimAr <- as.data.frame(TwoDimAr)
  names(TwoDimAr) <- c('ni', 'fit', 'fit.1', 'fit.2', 'fit.3')
  TwoDimAr$SR <- successRuns
  TwoDimAr$Nl1 <- Nl[1]; TwoDimAr$Nl2 <- Nl[2]; TwoDimAr$Nl3 <- Nl[3]
  TwoDimAr$Ng1 <- Ng[1]; TwoDimAr$Ng2 <- Ng[2]; TwoDimAr$Ng3 <- Ng[3]
  return(TwoDimAr)
}


flt <- function(scenario, Nl.vec, Ng.vec, nDay, gen.pd) {
  scenario <- cleanUpZeros(scenario)
  flatScenario <- flaten2D(scenario, Nl = Nl.vec, Ng = Ng.vec)
  flatScenario$Day <- rep(rep(1:nDay, each = gen.pd), flatScenario$SR[1])
  flatScenario.last <- flatScenario[seq(gen.pd, nrow(flatScenario), gen.pd), ] ## Keep only last generation per day
  
  days <- sort(unique(flatScenario.last$Day))
  fit <- tapply(flatScenario.last$fit, flatScenario.last$Day, mean)
  fit.1 <- tapply(flatScenario.last$fit.1, flatScenario.last$Day, mean)
  fit.2 <- tapply(flatScenario.last$fit.2, flatScenario.last$Day, mean)
  fit.3 <- tapply(flatScenario.last$fit.3, flatScenario.last$Day, mean)
  
  d <- cbind(fit, fit.1, fit.2, fit.3, days, 
             matrix(rep(c(Nl.vec, Ng.vec), nDay), nrow=nDay, ncol= length(c(Nl.vec, Ng.vec)), byrow=T))
  d <- as.data.frame(d)
  names(d)[6:11] <- c("Nl1", "Nl2", "Nl3", "Ng1", "Ng2", "Ng3")
  d
}


# All scenarios with same maxPsize = 2000 (modeling populations within day)
allScenarios.Day <- c("Nl.11015.Ng.101010",
                      "Nl.11015.Ng.21015", 
                      "Nl.11015.Ng.61012", 
                      "Nl.124.Ng.101010",  
                      "Nl.124.Ng.21015",   
                      "Nl.124.Ng.61012",   
                      "Nl.1610.Ng.101010", 
                      "Nl.1610.Ng.21015",  
                      "Nl.1610.Ng.61012")


 y1 <- as.data.frame(apply(Nl.11015.Ng.101010, c(1,2), mean))
 colnames(y1) <- c("ni", "fit", "fit.1", "fit.2", "fit.3")
 y1$Nl1 <- 1; y1$Nl2 <- 10; y1$Nl3 <- 15; y1$Ng1 <- 10; y1$Ng2 <- 10; y1$Ng3 <- 10; y1$SR <- 76
 y1$gen.number <- rep(1:24); y1$Day <- rep(1:20, each = 24)

# y2 <- as.data.frame(apply(Nl.11015.Ng.21015, c(1,2), mean))
# colnames(y2) <- c("ni", "fit", "fit.1", "fit.2", "fit.3")
# y2$Nl1 <- 1; y2$Nl2 <- 10; y2$Nl3 <- 15; y2$Ng1 <- 2; y2$Ng2 <- 10; y2$Ng3 <- 15; y2$SR <- 99
# y2$gen.number <- rep(1:24); y2$Day <- rep(1:20, each = 24)
# 
# y3 <- as.data.frame(apply(Nl.11015.Ng.61012, c(1,2), mean))
# colnames(y3) <- c("ni", "fit", "fit.1", "fit.2", "fit.3")
# y3$Nl1 <- 1; y3$Nl2 <- 10; y3$Nl3 <- 15; y3$Ng1 <- 6; y3$Ng2 <- 10; y3$Ng3 <- 12; y3$SR <- 61
# y3$gen.number <- rep(1:24); y3$Day <- rep(1:20, each = 24)
# 
# y4 <- as.data.frame(apply(Nl.124.Ng.101010, c(1,2), mean))
# colnames(y4) <- c("ni", "fit", "fit.1", "fit.2", "fit.3")
# y4$Nl1 <- 1; y4$Nl2 <- 2; y4$Nl3 <- 4; y4$Ng1 <- 10; y4$Ng2 <- 10; y4$Ng3 <- 10; y4$SR <- 75
# y4$gen.number <- rep(1:24); y4$Day <- rep(1:20, each = 24)
# 
# y5 <- as.data.frame(apply(Nl.124.Ng.21015, c(1,2), mean))
# colnames(y5) <- c("ni", "fit", "fit.1", "fit.2", "fit.3")
# y5$Nl1 <- 1; y5$Nl2 <- 2; y5$Nl3 <- 4; y5$Ng1 <- 2; y5$Ng2 <- 10; y5$Ng3 <- 15; y5$SR <- 99
# y5$gen.number <- rep(1:24); y5$Day <- rep(1:20, each = 24)
# 
# y6 <- as.data.frame(apply(Nl.124.Ng.61012, c(1,2), mean))
# colnames(y6) <- c("ni", "fit", "fit.1", "fit.2", "fit.3")
# y6$Nl1 <- 1; y6$Nl2 <- 2; y6$Nl3 <- 4; y6$Ng1 <- 6; y6$Ng2 <- 10; y6$Ng3 <- 12; y6$SR <- 60
# y6$gen.number <- rep(1:24); y6$Day <- rep(1:20, each = 24)
# 
# y7 <- as.data.frame(apply(Nl.1610.Ng.101010, c(1,2), mean))
# colnames(y7) <- c("ni", "fit", "fit.1", "fit.2", "fit.3")
# y7$Nl1 <- 1; y7$Nl2 <- 6; y7$Nl3 <- 10; y7$Ng1 <- 10; y7$Ng2 <- 10; y7$Ng3 <- 10; y7$SR <- 80
# y7$gen.number <- rep(1:24); y7$Day <- rep(1:20, each = 24)
# 
# y8 <- as.data.frame(apply(Nl.1610.Ng.21015, c(1,2), mean))
# colnames(y8) <- c("ni", "fit", "fit.1", "fit.2", "fit.3")
# y8$Nl1 <- 1; y8$Nl2 <- 6; y8$Nl3 <- 10; y8$Ng1 <- 2; y8$Ng2 <- 10; y8$Ng3 <- 15; y8$SR <- 100
# y8$gen.number <- rep(1:24); y8$Day <- rep(1:20, each = 24)
# 
# y9 <- as.data.frame(apply(Nl.1610.Ng.61012, c(1,2), mean))
# colnames(y9) <- c("ni", "fit", "fit.1", "fit.2", "fit.3")
# y9$Nl1 <- 1; y9$Nl2 <- 6; y9$Nl3 <- 10; y9$Ng1 <- 6; y9$Ng2 <- 10; y9$Ng3 <- 12; y9$SR <- 65
# y9$gen.number <- rep(1:24); y9$Day <- rep(1:20, each = 24)

# y10 <- as.data.frame(apply(Nl.11010.Ng.1055, c(1,2), mean))
# colnames(y10) <- c("ni", "fit", "fit.1", "fit.2", "fit.3")
# y10$Nl1 <- 1; y10$Nl2 <- 10; y10$Nl3 <- 10; y10$Ng1 <- 10; y10$Ng2 <- 5; y10$Ng3 <- 5; y10$SR <- 71
# y10$gen.number <- rep(1:24); y10$Day <- rep(1:20, each = 24)


y1 <- as.data.frame(flaten2D(Nl.11015.Ng.101010, Nl = c(1,10,15), Ng = c(10,10,10)))
y1$SR <- 76; y1$gen.number <- rep(1:24); y1$Day <- rep(1:20, each = 24)

y2 <- as.data.frame(flaten2D(Nl.11015.Ng.21015, Nl = c(1,10,15), Ng = c(2,10,15)))
y2$SR <- 99; y2$gen.number <- rep(1:24); y2$Day <- rep(1:20, each = 24)

y3 <- as.data.frame(flaten2D(Nl.11015.Ng.61012, Nl = c(1,10,15), Ng = c(6,10,12)))
y3$SR <- 61; y3$gen.number <- rep(1:24); y3$Day <- rep(1:20, each = 24)

y4 <- as.data.frame(flaten2D(Nl.124.Ng.101010, Nl = c(1,2,4), Ng = c(10,10,10)))
y4$SR <- 75; y4$gen.number <- rep(1:24); y4$Day <- rep(1:20, each = 24)

y5 <- as.data.frame(flaten2D(Nl.124.Ng.21015, Nl = c(1,2,4), Ng = c(2,10,15)))
y5$SR <- 99; y5$gen.number <- rep(1:24); y5$Day <- rep(1:20, each = 24)

y6 <- as.data.frame(flaten2D(Nl.124.Ng.61012, Nl = c(1,2,4), Ng = c(6,10,12)))
y6$SR <- 60; y6$gen.number <- rep(1:24); y6$Day <- rep(1:20, each = 24)

y7 <- as.data.frame(flaten2D(Nl.1610.Ng.101010, Nl = c(1,6,10), Ng = c(10,10,10)))
y7$SR <- 80; y7$gen.number <- rep(1:24); y7$Day <- rep(1:20, each = 24)

y8 <- as.data.frame(flaten2D(Nl.1610.Ng.21015, Nl = c(1,6,10), Ng = c(2,10,15)))
y8$SR <- 100; y8$gen.number <- rep(1:24); y8$Day <- rep(1:20, each = 24)

y9 <- as.data.frame(flaten2D(Nl.1610.Ng.61012, Nl = c(1,6,10), Ng = c(6,10,12)))
y9$SR <- 65; y9$gen.number <- rep(1:24); y9$Day <- rep(1:20, each = 24)


genDay <- rbind(y1,y2,y3,y4,y5,y6,y7,y8,y9)

fitGen.lm <- lm(fit ~ poly(gen.number,2, raw = T) + Nl2 + Nl3 + Ng1 + Ng3, genDay[genDay$Day == 2,]) #input whatever day trying to predict 
dataGenFit <- data.frame(Day = 2, Nl2 = 6, Nl3 = 10, Ng1 = 10, Ng3 = 10, gen.number = 1:24)
fitGen.pred <- predict(fitGen.lm, newdata = dataGenFit)

genDay.lm <- lm(ni ~ poly(gen.number,3, raw = T) + Nl2 + Nl3 + Ng1 + Ng3 + fit, genDay[genDay$Day == 2,]) #input whatever day trying to predict 
dataGenDay <- data.frame(Day = 2, Nl2 = 6, Nl3 = 10, Ng1 = 10, Ng3 = 10, gen.number = 1:24, fit = fitGen.pred)
genDay.pred <- predict(genDay.lm, newdata = dataGenDay)





#plotting y_ generations/day results
plotSigmoid <- function(avgScenario, gen.pd, nDay, tt){
  plot(avgScenario$gen.number, avgScenario$ni, main = tt)
  lines(avgScenario[1:24,]$gen.number, predict(lm(ni ~ poly(gen.number, 3), avgScenario[avgScenario$Day == 1,])))
  for(i in 1:nDay){
    lines(avgScenario[(i*24+1):((i+1)*24),]$gen.number, avgScenario[(i*24+1):((i+1)*24),]$ni)
  }
}

plotSigmoid(y1,24,20, "Nl.11015.Ng.101010")
plotSigmoid(y2,24,20, "Nl.11015.Ng.21015")
plotSigmoid(y3,24,20, "Nl.11015.Ng.61012")
plotSigmoid(y4,24,20, "Nl.124.Ng.101010")
plotSigmoid(y5,24, 20, "Nl.124.Ng.21015")
plotSigmoid(y6,24, 20, "Nl.124.Ng.61012")
plotSigmoid(y7,24, 20, "Nl.1610.Ng.101010")
plotSigmoid(y8,24, 20, "Nl.1610.Ng.21015")
plotSigmoid(y9,24, 20, "Nl.1610.Ng.61012")


##Parent sigmoid
sigmoid <- function(x){
  100 / (1 + exp(-x+12))
}
plot(1:24, sigmoid(1:24), type = "l", lty = 2, lwd = 3, col = "red")



allScenarios <- c("Nl.11015.Ng.101010",
                  "Nl.11015.Ng.21015", 
                  "Nl.11015.Ng.61012", 
                  "Nl.124.Ng.101010",  
                  "Nl.124.Ng.21015",   
                  "Nl.124.Ng.61012",   
                  "Nl.1610.Ng.101010", 
                  "Nl.1610.Ng.21015",  
                  "Nl.1610.Ng.61012",
                  "Nl.11010.Ng.1055",  #Rm = 0.003, maxPsize = 10000
                  "Nl.11515.Ng.1055",  #Rm = 0.003, maxPsize = 10000
                  #"Nl.11515.Ng.1055.Rm5",
                  "Nl.151515.Ng.1055") #Rm = 0.003, maxPsize = 10000

x1 <- flt(Nl.124.Ng.21015, c(1,2,4), c(2,10,15), nDay = 20, gen.pd = 24)
x2 <- flt(Nl.124.Ng.61012, c(1,2,4), c(6,10,12), nDay = 20, gen.pd = 24)
x3 <- flt(Nl.124.Ng.101010, c(1,2,4), c(10,10,10), nDay = 20, gen.pd = 24)

x4 <- flt(Nl.1610.Ng.21015, c(1,6,10), c(2,10,15), nDay = 20, gen.pd = 24)
x5 <- flt(Nl.1610.Ng.61012, c(1,6,10), c(6,10,12), nDay = 20, gen.pd = 24)
x6 <- flt(Nl.1610.Ng.101010, c(1,6,10), c(10,10,10), nDay = 20, gen.pd = 24)

x7 <- flt(Nl.11015.Ng.21015, c(1,10,15), c(2,10,15), nDay = 20, gen.pd = 24)
x8 <- flt(Nl.11015.Ng.61012, c(1,10,15), c(6,10,12), nDay = 20, gen.pd = 24)
x9 <- flt(Nl.11015.Ng.101010, c(1,10,15), c(10,10,10), nDay = 20, gen.pd = 24)

x10 <- flt(Nl.11010.Ng.1055, c(1,10,10), c(10,5,5), nDay = 20, gen.pd = 24)
x11 <- flt(Nl.11515.Ng.1055, c(1,15,15), c(10,5,5), nDay = 20, gen.pd = 24) 
x12 <- flt(Nl.151515.Ng.1055, c(15,15,15), c(10,5,5), nDay = 20, gen.pd = 24)


d <- rbind(x1,x2,x3,x4,x5,x6,x7,x8,x9, x10, x11, x12)
# Construct d with an end field of scenario number
  # d$scenario.ms <- rep(1:3, each = 60)
  # d$scenario.g <- rep(rep(1:3, each = 20), 3)

# Construct d with an end field of success rate
d$success <- 0
d[d$Nl1 == 1 & d$Nl2 == 2 & d$Nl3 == 4 & d$Ng1 == 2 & d$Ng2 == 10 & d$Ng3 == 15,"success"] <- 99
d[d$Nl1 == 1 & d$Nl2 == 2 & d$Nl3 == 4 & d$Ng1 == 6 & d$Ng2 == 10 & d$Ng3 == 12,"success"] <- 60
d[d$Nl1 == 1 & d$Nl2 == 2 & d$Nl3 == 4 & d$Ng1 == 10 & d$Ng2 == 10 & d$Ng3 == 10,"success"] <- 75
d[d$Nl1 == 1 & d$Nl2 == 6 & d$Nl3 == 10 & d$Ng1 == 2 & d$Ng2 == 10 & d$Ng3 == 15,"success"] <- 100
d[d$Nl1 == 1 & d$Nl2 == 6 & d$Nl3 == 10 & d$Ng1 == 6 & d$Ng2 == 10 & d$Ng3 == 12,"success"] <- 65
d[d$Nl1 == 1 & d$Nl2 == 6 & d$Nl3 == 10 & d$Ng1 == 10 & d$Ng2 == 10 & d$Ng3 == 10,"success"] <- 80
d[d$Nl1 == 1 & d$Nl2 == 10 & d$Nl3 == 15 & d$Ng1 == 2 & d$Ng2 == 10 & d$Ng3 == 15,"success"] <- 99
d[d$Nl1 == 1 & d$Nl2 == 10 & d$Nl3 == 15 & d$Ng1 == 6 & d$Ng2 == 10 & d$Ng3 == 12,"success"] <- 61
d[d$Nl1 == 1 & d$Nl2 == 10 & d$Nl3 == 15 & d$Ng1 == 10 & d$Ng2 == 10 & d$Ng3 == 10,"success"] <- 76
d[d$Nl1 == 1 & d$Nl2 == 10 & d$Nl3 == 10 & d$Ng1 == 10 & d$Ng2 == 5 & d$Ng3 == 5,"success"] <- 71
d[d$Nl1 == 1 & d$Nl2 == 15 & d$Nl3 == 15 & d$Ng1 == 10 & d$Ng2 == 5 & d$Ng3 == 5,"success"] <- 63
d[d$Nl1 == 15 & d$Nl2 == 15 & d$Nl3 == 15 & d$Ng1 == 10 & d$Ng2 == 5 & d$Ng3 == 5,"success"] <- 11


d.lm1 <- lm(fit.1 ~ poly(days,3, raw = T) + Nl1 + Nl2 + Nl3 + Ng1 + Ng2 + Ng3, d)
d.lm2 <- lm(fit.2 ~ days + Nl1 + Nl2 + Nl3 + Ng1 + Ng2 + Ng3, d)
d.lm3 <- lm(fit.3 ~ poly(days,3, raw = T) + Nl1 + Nl2 + Nl3 + Ng1 + Ng2 + Ng3, d)

# Introduce scenario and interaction to the LM
d.lm1 <- lm(fit.1 ~ poly(days,3, raw = T) + Nl2 + Nl3 + Ng1 + Ng3 + scenario.ms:scenario.g, d)
d.lm2 <- lm(fit.2 ~ poly(days,2, raw = T) + Nl2 + Nl3 + Ng1 + Ng3 + scenario.ms:scenario.g, d)
d.lm3 <- lm(fit.3 ~ poly(days,2, raw = T) + Nl2 + Nl3 + Ng1 + Ng3 + scenario.ms:scenario.g, d)

# Predict success rate from sections
success.lm <- lm(success ~ Nl1 + Nl2 + Nl3 + Ng1 + Ng2 + Ng3, d)

# User input
Nl1 <- 3; Nl2 <- 6; Nl3 = 10; Ng1 = 6; Ng2 = 10; Ng3 = 12; scenario.ms = 2; scenario.g = 2

# construct new data using user input
newD <- data.frame(days = 1:30, Nl2, Nl3, Ng1, Ng3, scenario.ms, scenario.g)
newD1 <- data.frame(days = 1:30, Nl2, Nl3, Ng1, Ng3, scenario.ms, scenario.g)
newD2 <- data.frame(days = 1:30, Nl2, Nl3, Ng1, Ng3, scenario.ms, scenario.g)
newD3 <- data.frame(days = 1:30, Nl2, Nl3, Ng1, Ng3, scenario.ms, scenario.g)

# predict fitness based on user input 
fit1.pred <- predict(d.lm1, newdata = newD1)
for(i in 1:length(fit1.pred)){
  if (fit1.pred[i] > 1) fit1.pred[i] = 1
}

fit2.pred <- predict(d.lm2, newdata = newD2)
for(i in 1:length(fit2.pred)){
  if (fit2.pred2[i] > 1) fit2.pred[i] = 1
}

fit3.pred <- predict(d.lm3, newdata = newD3)
for(i in 1:length(fit3.pred)){
  if (fit3.pred[i] > 1) fit3.pred[i] = 1
}

fit.pred <- (fit1.pred + fit2.pred + fit3.pred)/3


plot(1:30, fit1.pred)
points(1:30, fit2.pred)
points(1:30, fit3.pred)
points(1:30, fit.pred)







