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

allScenarios <- c("Nl.11015.Ng.101010",
                  "Nl.11015.Ng.21015", 
                  "Nl.11015.Ng.61012", 
                  "Nl.124.Ng.101010",  
                  "Nl.124.Ng.21015",   
                  "Nl.124.Ng.61012",   
                  "Nl.1610.Ng.101010", 
                  "Nl.1610.Ng.21015",  
                  "Nl.1610.Ng.61012",
                  "Nl.3612.Ng.21015")

x1 <- flt(Nl.124.Ng.21015, c(1,2,4), c(2,10,15), nDay = 20, gen.pd = 24)
x2 <- flt(Nl.124.Ng.61012, c(1,2,4), c(6,10,12), nDay = 20, gen.pd = 24)
x3 <- flt(Nl.124.Ng.101010, c(1,2,4), c(10,10,10), nDay = 20, gen.pd = 24)

x4 <- flt(Nl.1610.Ng.21015, c(1,6,10), c(2,10,15), nDay = 20, gen.pd = 24)
x5 <- flt(Nl.1610.Ng.61012, c(1,6,10), c(6,10,12), nDay = 20, gen.pd = 24)
x6 <- flt(Nl.1610.Ng.101010, c(1,6,10), c(10,10,10), nDay = 20, gen.pd = 24)

x7 <- flt(Nl.11015.Ng.21015, c(1,10,15), c(2,10,15), nDay = 20, gen.pd = 24)
x8 <- flt(Nl.11015.Ng.61012, c(1,10,15), c(6,10,12), nDay = 20, gen.pd = 24)
x9 <- flt(Nl.11015.Ng.101010, c(1,10,15), c(10,10,10), nDay = 20, gen.pd = 24)

# x10 <- flt(Nl.3612.Ng.21015, c(3,6,12), c(3,6,12), nDay = 20, gen.pd = 24)
# x13 <- flt(Nl.124.Ng.4816, c(1,2,4), c(4,8,16), nDay = 20, gen.pd = 24)

d <- rbind(x1,x2,x3,x4,x5,x6,x7,x8,x9)
# Construct d with an end field of scenario number
d$scenario.ms <- rep(1:3, each = 60)
d$scenario.g <- rep(rep(1:3, each = 20), 3)
# Construct d with an end field of success rate
d$success <- 0
d[d$scenario.ms == 1 & d$scenario.g == 1,"success"] <- 99
d[d$scenario.ms == 1 & d$scenario.g == 2,"success"] <- 60
d[d$scenario.ms == 1 & d$scenario.g == 3,"success"] <- 75
d[d$scenario.ms == 2 & d$scenario.g == 1,"success"] <- 100
d[d$scenario.ms == 2 & d$scenario.g == 2,"success"] <- 65
d[d$scenario.ms == 2 & d$scenario.g == 3,"success"] <- 80
d[d$scenario.ms == 3 & d$scenario.g == 1,"success"] <- 99
d[d$scenario.ms == 3 & d$scenario.g == 2,"success"] <- 61
d[d$scenario.ms == 3 & d$scenario.g == 3,"success"] <- 76



d.lm1 <- lm(fit.1 ~ poly(days,2, raw = T) + Nl2 + Nl3 + Ng1 + Ng3, d)
d.lm2 <- lm(fit.2 ~ poly(days,2, raw = T) + Nl2 + Nl3 + Ng1 + Ng3, d)
d.lm3 <- lm(fit.3 ~ poly(days,2, raw = T) + Nl2 + Nl3 + Ng1 + Ng3, d)

# Introduce scenario and interaction to the LM
d.lm1 <- lm(fit.1 ~ poly(days,3, raw = T) + Nl2 + Nl3 + Ng1 + Ng3 + scenario.ms:scenario.g, d)
d.lm2 <- lm(fit.2 ~ poly(days,2, raw = T) + Nl2 + Nl3 + Ng1 + Ng3 + scenario.ms:scenario.g, d)
d.lm3 <- lm(fit.3 ~ poly(days,2, raw = T) + Nl2 + Nl3 + Ng1 + Ng3 + scenario.ms:scenario.g, d)

# Predict success rate from sections
success.lm <- lm(success ~ Nl2 + Nl3 + Ng1 + Ng3, d)

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







