### TSIR for small populations

library(plotrix)
library(Rwave)
library(MASS)
library(signal)
library(stats)
library(polynom)
library(quantchem)
library(matrixcalc)
library(stats)
library(e1071)

place <- c("reykjavik", "hafnafjordur", "akureyri","bornholm", "faroe")
folder <- c("iceland", "iceland","iceland","bornholm", "faroe")
name <- 1

# Set up parameters
periodicity <- 24
vaccine <- 1965
delay <- 8
sensitivity <- 5
sensitivity.length <- 3
numsim <- 1000

# Import data
setwd("~/Documents/Grenfell Research/Measles")
data <- read.csv(paste("SmallPopTSIR/data/", folder[name], "/",place[name], ".csv", sep = ""))
#data <- read.csv("SmallPopTSIR/data/faroe/faroe.csv")
#data <- read.csv("SmallPopTSIR/data/bornholm/bornholm.csv")
plot(data$reported_cases, type = "l")
data <- data[data$time <= vaccine, ]


time <- data$time
lengthdata <- length(time)

B <- data$births
C <- data$reported_cases

png(file=paste("~/Documents/Grenfell Research/Measles/SmallPopTSIR/data/", folder[name], "/results/", place[name],"_0_timeseries_cases.png", sep = ""),
    width=800,height=700)
plot(time, C, type = "l")
dev.off()

png(file=paste("~/Documents/Grenfell Research/Measles/SmallPopTSIR/data/", folder[name], "/results/", place[name],"_0_timeseries_births.png", sep = ""),
    width=800,height=700)
plot(time, B, type = "l")
dev.off()


## Notation from paper
X <- cumsum(C) # cumulative cases
Y <- cumsum(B) # cumulative births

plot(Y, type = "l", ylab = "Cumulative Births", xlab = "Time")
plot(X, type = "l", ylab = "Cumulative Cases", xlab = "Time")

##################### Susceptible reconstruction ###########################################

### Lowess regression
lowess.fit <- lowess(X, Y, f = 2/3, iter = 0)
plot(time, Y, type = "l", ylab = "Cumulative Births/Cases")
lines(time, lowess.fit$y, type = "l", col = "blue") #predicted cumulative births (blue)
lines(time, cumsum(C), col = "red")

plot(X, Y, type = "l")
lines(lowess.fit$x, lowess.fit$y, type = "l")

spline <- smooth.spline(lowess.fit$x,  lowess.fit$y, df = 2.5)

Yhat <- lowess.fit$y
rho <- (predict(spline, deriv=1)$y)
Z <- -(Yhat - Y)


############################

# polynomial
cum.reg <- lm(Y ~ X + I(X^2) + I(X^3)) 

Yhat <- predict(cum.reg)
rho <- derivative(cum.reg, X)
m.rho <- round(mean(1/rho), 2)
Z <- -(Yhat - Y)

png(file=paste("~/Documents/Grenfell Research/Measles/SmallPopTSIR/data/", folder[name], "/results/", place[name],"_1_Yhat.png", sep = ""),
    width=800,height=700)
plot(time, Y, type = "l", ylab = "Cumulative Births/Cases", col = "green",
     main = "Cumulative Births (green), Inferred Cumulative Births (red), 
     \n Reported Cumulative Cases (blue)",
     cex.main = 1)
lines(time, Yhat, type = "l", col = "red") #predicted cumulative births (red)
lines(time, cumsum(C), col = "blue")
dev.off()


png(file=paste("~/Documents/Grenfell Research/Measles/SmallPopTSIR/data/", folder[name], "/results/", place[name],"_1_rho.png", sep = "")
    ,width=800,height=700)
plot(time, 1/rho, type = "l", 
     main = paste("Inferred Reporting Rate, mean 1/rho = ", m.rho, sep = ""), cex.main = 1)
dev.off()

png(file=paste("~/Documents/Grenfell Research/Measles/SmallPopTSIR/data/", folder[name], "/results/", place[name],"_1_Z.png", sep = "")
    ,width=800,height=700)
plot(time, Z , type = "l", xlab = "Time", ylab = "Residual", main = "Susceptible Dynamics, Z(t)", 
     cex.main = 1) #plot Z
dev.off()

# Ic are actual cases: reported cases multiplied by rho
Ic <- C*rho 
plot(Ic, type = "l")
lines(C, type = "l", col = "red")

## Keep only positive values
keep.index <- which(Ic[1:lengthdata - 1] > sensitivity & Ic[2:lengthdata] > sensitivity)
length.pos <- length(keep.index)

## Set up the parameters
seas = rep(1:periodicity, lengthdata/periodicity)[1:lengthdata - 1][keep.index]
lInew <- log(Ic[2:lengthdata]) [keep.index]
lIold <- log(Ic[1:lengthdata -1]) [keep.index]
Zold = Z[1:lengthdata - 1] [keep.index]


######### This is from Bjornstadt (2007): Susceptible Reconstruction #######################
#Smean = seq(-0.015, 0.2, by=0.001)*.2E5  
min(Zold)
Smean = seq(abs(min(Zold)) + 1, abs(min(Zold))*4, by = 100) ## Reasonable range of candidate values
llik = rep(NA, length(Smean))
summary(Smean)

## Loop over all candidate values of Smean. Glm is used with no intercept so that as.factor(seas) becomes the 
## estimates of the log-beta's. glmfit$deviance holds -2*log-likelihood

for (i in 1:length(Smean)) {
  lSold = log(Smean[i] + Zold)
  glmfit = glm(lInew ~ -1 + as.factor(seas) + lIold + offset(lSold)) 
  # -1 removes intercept, offset = known coeff of 1
  #llik[i] = glmfit$deviance # -2*loglikelihood
  llik[i] = logLik(glmfit)
}
Sbar <- Smean[which(llik == max(llik))] #we want to maximize the log likelihood
Sbar

png(file=paste("~/Documents/Grenfell Research/Measles/SmallPopTSIR/data/", folder[name], "/results/", place[name],"_2_Sbar.png", sep = "")
    ,width=800,height=700)
plot(Smean, -llik, type = "l", log = "y", main = paste("Sbar = ", round(Sbar, 2), sep = ""))
abline(v=Sbar)
dev.off()

lSold = log(Sbar + Zold) ## This is the best estimate (maximizes log likelihood)
# best estimates
lmfit <- lm(lInew ~ -1 + as.factor(seas) + lIold + offset(lSold))
sum <- summary(lmfit)
alpha <- lmfit$coef[periodicity + 1]


######################################################

#### Check other method: Taylor expansion
lmfit2 <- lm(lInew ~ -1 + as.factor(seas) + lIold + Zold)
sum2 <- summary(lmfit2)

png(file=paste("~/Documents/Grenfell Research/Measles/SmallPopTSIR/data/", folder[name], "/results/", place[name],"_2_Seasonality_compare.png", sep = "")
    ,width=800,height=700)
plot(1:periodicity, exp(lmfit2$coef[1:periodicity]), type = "l", col = "red", 
     ylab = "Seasonal Coeff", xlab = "Week Number", main = "Seasonality")
lines(1:periodicity, exp(lmfit$coef[1:periodicity]) * Sbar, type = "l", lty = 2, col = "blue")
legend("bottomleft", legend = c("Using ML", "Using Taylor Exp"), col = c("blue", "red"), lty = c(2,1), cex = 0.5)
dev.off()

# Estimates from this method (lmfit2)
alpha2 <- lmfit2$coef[periodicity + 1] #alfa
zeta <- lmfit2$coef[periodicity + 2] #zeta
mean(lmfit2$coef[1:periodicity])
Sbar2 <- 1/zeta
#####################

lrstar <- rep(lmfit$coef[1:periodicity], lengthdata/periodicity)
length(lrstar)


# Confidence intervals
se <- sum$coefficients[ ,2]
upper <- lmfit$coef[1:periodicity] - mean(lmfit$coef[1:periodicity]) + qnorm(0.975) * se[1:periodicity]
lower <- lmfit$coef[1:periodicity] - mean(lmfit$coef[1:periodicity]) - qnorm(0.975) * se[1:periodicity]

# Plotting mean-centered seasonality
png(file=paste("~/Documents/Grenfell Research/Measles/SmallPopTSIR/data/", folder[name], "/results/", place[name],"_2_Seasonality.png", sep = "")
    ,width=800,height=700)
plot(1:periodicity, lmfit$coef[1:periodicity] - mean(lmfit$coef[1:periodicity]), 
     type = "l", ylab = "Seasonality", ylim = c(min(lower), max(upper)), 
     xaxp = c(0, periodicity+2, (periodicity+2)/2), xlab = "Period")
abline(h = 0)
lines(1:periodicity, upper, type = "l", lty = 2)
lines(1:periodicity, lower, type = "l", lty = 2)
dev.off()


################## Simulation #######################################


### Identify epidemics 
#look at sequences above sensitivity
start.index <- ifelse(keep.index[2:length.pos] - keep.index[1:length.pos - 1] > 1, keep.index[2:length.pos], NA)
start.index <- append(keep.index[1], start.index)
start.index <- start.index[is.na(start.index) == FALSE]
end.index <- ifelse(keep.index[2:length.pos] - keep.index[1:length.pos - 1] > 1, keep.index[1:length.pos - 1], NA)
end.index <- end.index[is.na(end.index) == FALSE]
end.index <- append(end.index, keep.index[length(keep.index)])

#merge two sequences if the end of last sequence and beginning of next sequence differ by less than 3
diff <- which(start.index[2:length(start.index)] - end.index[1:length(end.index) - 1] <= sensitivity.length)
for (i in 1:length(diff)){
  j <- diff[i]
  start.index[j + 1] <- NA
  end.index[j] <- ifelse(is.na(end.index[j]) == FALSE, end.index[j + 1], NA)
  end.index[j + 1] <- NA
}
start.index <- start.index[is.na(start.index) == FALSE]
end.index <- end.index[is.na(end.index) == FALSE]


plot(Ic, type = "l")
for (i in 1:length(start.index)) {
  abline(v = end.index[i], col = "blue")
  abline(v = start.index[i], col = "green")
}



# real duration 
duration <- end.index - start.index


#create a matrix for each epidemic with start index, end index and duration
epi <- as.data.frame(as.matrix(cbind(start.index, end.index, duration)))
epi <- epi[is.na(epi$duration) == FALSE, ]
epi <- epi[epi$duration > sensitivity.length, ] #only keep epidemics that are greater than a specified length

# real size
for (e in 1:length(epi$start.index)) {
  epi$size[e] <- sum(Ic[epi$start.index[e]: epi$end.index[e]])
}


# logI and logS are the full time series vectors
logI <- rep(0, lengthdata)
logI[keep.index] <- log(Ic[keep.index])
logS <- rep(0, lengthdata)
logS[keep.index] = log(Sbar + Z[keep.index])

#S0
epi$S0 <- exp(logS[epi$start.index])


# Set up parameters
predS <- matrix(NA, lengthdata, numsim)
predI <- matrix(NA, lengthdata, numsim)
lambda <- matrix(NA, lengthdata, numsim)
rstar <- exp(lrstar)
alpha <- 0.97


for (n in 1:numsim) {
  
  for (i in 1: length(epi$start.index)) {
    start <- epi$start.index[i]
    end <- ifelse( i < length(epi$start.index), epi$start.index[i+1] - 1, lengthdata)
    
    # initial conditions
    predS[start, n] <- exp(logS[start]) #This is where our inferred Sbar comes in
    predI[start, n] <- exp(logI[start])
    
    
    for (j in 1:(end - start)) {
      t <- start + j
      
      
      lambda[t, n] <- rstar[t] * (predI[(t-1), n]^alpha) * predS[(t-1), n]
      if(lambda[t, n] > 0 & t != lengthdata) {
        predI[t, n] <- rnegbin(1, lambda[t,n], max(round(predI[(t - 1), n]), 1))
      }
      else {
        predI[t, n] <- 0
      }
      
      predS[t,n] <- round(B[t - delay] + predS[t-1, n] - predI[t, n]) 
    }
    
  }
  
}

# Transpose the matrix so that rows = simulations, and columns = time
predI <- t(predI)

I <- colMeans(predI)
S <- colMeans(predS)
I <- ifelse(is.na(I) == TRUE, 0, I)
plot(I, type = "l")

png(file=paste("~/Documents/Grenfell Research/Measles/SmallPopTSIR/data/", folder[name], "/results/", place[name],"_3_predictions_R.png", sep = "")
    ,width=800,height=700)
plot(time, predI[1,], type = "l", col = "pink", ylim = c(0, max(Ic) + 1500), ylab = "Cases",
     main = paste("Predicted cases from ", numsim, " simulations (red) 
        Actual cases corrected for underreporting (blue)", sep = ""), 
     cex.main = 0.7)
for (i in 2:nrow(predI)) {
  lines(time, predI[i, ], type = "l", col = "pink")
}
lines(time, I, type = "l", col = "red", lwd = 2)
lines(time, Ic, type = "l", col = "blue")
dev.off()

## Non-parametric bootstrapping for confidence intervals for prediction
numboot <- 1000
matrix.boot <- matrix(NA, ncol = lengthdata, nrow = numboot)


for (i in 1:numboot) {
  boot.smple <- predI[sample(numsim, replace = TRUE), ] #Take 1000 simulations with replacement
  matrix.boot[i,] <- colMeans(boot.smple) #calculate mean
}


lower.bound <- apply(matrix.boot, 2, function(x) quantile(x, 0.025, na.rm = TRUE))
upper.bound <- apply(matrix.boot, 2, function(x) quantile(x, 0.975, na.rm = TRUE))


png(file=paste("~/Documents/Grenfell Research/Measles/SmallPopTSIR/data/", folder[name], "/results/", place[name],"_3_predictions_CI_R.png", sep = "")
    ,width=800,height=700)
plot(time, I, type = "l", col = "red", ylim = c(0, max(Ic) + 100), ylab = "Cases",
     main = paste("Predicted cases from ", numsim, " simulations (red) 
        Actual cases corrected for underreporting (blue)", sep = ""), 
     cex.main = 0.7)
lines(time, lower.bound, type = "l", lty = 2)
lines(time, upper.bound, type = "l", lty = 2)
lines(time, Ic, type = "l", col = "blue")
dev.off()


# Make sure to transpose the matrix again; fix code below later in order to avoid this
predI <- t(predI)

### Real vs predicted duration and size
pred.dur <- matrix(NA, length(epi$start.index), numsim)
pred.endindex <- matrix(NA, length(epi$start.index), numsim)
pred.size <- matrix(NA, length(epi$start.index), numsim)

for (c in 1:ncol(predI)) {
  predIcol <- predI[,c]
  keep.index.sim <- which(predIcol[1:lengthdata - 1] > 0 & predIcol[2:lengthdata] > 0)
  length.pos.sim <- length(keep.index.sim)
  start.index.sim <- epi$start.index
  
  for (s in 1:length(start.index.sim)){
    zero <- which(predIcol[start.index.sim[s]:lengthdata] == 0)
    if (s != length(start.index.sim)) {
      pred.endindex[s,c] <- ifelse(min(zero) + start.index.sim[s] < start.index.sim[s+1], 
                                   (min(zero) + start.index.sim[s]) - 1, start.index.sim[s+1] - 1)
    } else {
      pred.endindex[s,c] <- (min(zero) + start.index.sim[s]) - 1
    }
    
    pred.dur[s,c] <- pred.endindex[s,c] - start.index.sim[s]
    pred.size[s,c] <- sum(predIcol[epi$start.index[s]: pred.endindex[s,c]])
    
  }
  
  
}


pred.dur <- as.data.frame(pred.dur)
pred.dur$real.dur <- epi$duration
pred.dur$id <- c(1:length(epi$start.index))
pred.dur$realS0 <- epi$S0

duration <- reshape(pred.dur, 
               varying = list(names(pred.dur)[1:numsim]), 
               v.names = "pred.dur",
               idvar = c("id", "real.dur", "realS0") ,     
               direction = "long")


# Real duration versus predicted duration
fit.dur <- lm(duration$pred.dur ~ -1 + duration$real.dur)
slope.dur <- round(fit.dur$coeff[1], 5)
rsq.dur <- round(summary(fit.dur)$r.squared, 5)
ci.dur <- predict(fit.dur, interval = "confidence", level = 0.95, newdata = duration)

png(file=paste("~/Documents/Grenfell Research/Measles/SmallPopTSIR/data/", folder[name], "/results/", place[name],"_4_Duration.png", sep = "")
    ,width=800,height=700)
plot(duration$real.dur, duration$pred.dur, pch = 1, cex = 0.5, 
     xlab = "Real Duration", ylab = "Predicted Duration", 
     main = paste("Slope = ", slope.dur, " R-squared = ", rsq.dur, sep = ""), cex.main = 0.9)
abline(fit.dur, col = "blue")
lines(duration$real.dur, ci.dur[ ,2], col = "red", lty = 3)
lines(duration$real.dur, ci.dur[ ,3], col = "red", lty = 3)
dev.off()


# S0 versus real duration
fitS0 <- lm(epi$S0 ~ epi$duration)
slope.S0.realdur <- round(fitS0$coeff[2], 5)
ci.S0.realdur <- predict(fitS0, interval = "confidence", level = 0.95, newdata = epi)


png(file=paste("~/Documents/Grenfell Research/Measles/SmallPopTSIR/data/", folder[name], "/results/", place[name],"_5_S0_real_dur.png", sep = "")
    ,width=800,height=700)
plot(epi$duration, epi$S0, pch = 1, cex = 0.5,
     xlab = "Real Duration", ylab = "S0", 
     main = paste("Slope = ", slope.S0.realdur, sep = ""), cex.main = 0.9)
abline(fitS0, col = "blue")
dev.off()

#S0 versus predicted size
fitS0.pred <- lm(duration$realS0 ~ duration$pred.dur)
slope.S0.preddur <- round(fitS0.pred$coeff[2], 5)
ci.S0.preddur <- predict(fitS0.pred, interval = "confidence", level = 0.95, newdata = duration)


png(file=paste("~/Documents/Grenfell Research/Measles/SmallPopTSIR/data/", folder[name], "/results/", place[name],"_5_S0_pred_dur.png", sep = "")
    ,width=800,height=700)
plot(duration$pred.dur, duration$realS0, pch = 1, cex = 0.5,
     xlab = "Predicted Duration", ylab = "S0", 
     main = paste("Slope = ", slope.S0.preddur, sep = ""), cex.main = 0.9)
abline(fitS0.pred, col = "blue")
dev.off()


### Real vs predicted size of each epidemic
pred.size <- as.data.frame(pred.size)
pred.size$real.size <- epi$size
pred.size$id <- c(1:length(epi$start.index))
pred.size$realS0 <- epi$S0

size <- reshape(pred.size, 
                varying = list(names(pred.size)[1:numsim]), 
                v.names = "pred.size",
                idvar = c("real.size", "id", "realS0"),      
                direction = "long")

# Real size versus predicted size
fit.size <- lm(size$pred.size ~ -1 + size$real.size)
summary(fit.size)
slope.size <- round(fit.size$coeff[1], 5)
rsq.size <- round(summary(fit.size)$r.squared, 5)
ci.size <- predict(fit.size, interval = "confidence", level = 0.95, newdata = size)


png(file=paste("~/Documents/Grenfell Research/Measles/SmallPopTSIR/data/", folder[name], "/results/", place[name],"_4_Size.png", sep = "")
    ,width=800,height=700)
plot(size$real.size, size$pred.size, pch = 1, cex = 0.5,
     xlab = "Real Size", ylab = "Predicted Size", 
     main = paste("Slope = ", slope.size, " R-squared = ", rsq.size, sep = ""), cex.main = 0.9)
abline(fit.size, col = "blue")
lines(size$real.size, ci.size[ ,2], col = "red", lty = 3)
lines(size$real.size, ci.size[ ,3], col = "red", lty = 3)
dev.off()

# S0 versus real size
fitS0 <- lm(epi$S0 ~ epi$size)
summary(fitS0)
slope.S0.realsize <- round(fitS0$coeff[2], 5)

png(file=paste("~/Documents/Grenfell Research/Measles/SmallPopTSIR/data/", folder[name], "/results/", place[name],"_5_S0_real.png", sep = "")
    ,width=800,height=700)
plot(epi$size, epi$S0, pch = 1, cex = 0.5,
     xlab = "Real Size", ylab = "S0", 
     main = paste("Slope = ", slope.S0.realsize, sep = ""), cex.main = 0.9)
abline(fitS0, col = "blue")
dev.off()

#S0 versus predicted size
fitS0.pred <- lm(size$realS0 ~ size$pred.size)
summary(fitS0.pred)
slope.S0.predsize <- round(fitS0.pred$coeff[2], 5)
ci.S0.pred <- predict(fitS0.pred, interval = "confidence", level = 0.95, newdata = size)

png(file=paste("~/Documents/Grenfell Research/Measles/SmallPopTSIR/data/", folder[name], "/results/", place[name],"_5_S0_pred.png", sep = "")
    ,width=800,height=700)
plot(size$pred.size, size$realS0, pch = 1, cex = 0.5,
     xlab = "Predicted Size", ylab = "S0", 
     main = paste("Slope = ", slope.S0.predsize, sep = ""), cex.main = 0.9)
abline(fitS0.pred, col = "blue")
lines(size$pred.size, ci.S0.pred[ ,2], col = "red", lty = 3)
lines(size$pred.size, ci.S0.pred[ ,3], col = "red", lty = 3)
dev.off()

