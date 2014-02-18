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

# Set up parameters
periodicity <- 24
vaccine <- 1965
delay <- 8
sensitivity <- 0
numsim <- 1000

# Import data
setwd("~/Documents/Grenfell Research/Measles")
data <- read.csv("SmallPopTSIR/data/iceland/reykjavik.csv")
#data <- read.csv("SmallPopTSIR/data/faroe/faroe.csv")
plot(data$reported_cases, type = "l")
data <- data[data$time <= vaccine, ]


time <- data$time
lengthdata <- length(time)

B <- data$births
plot(time, B, type = "l")


C <- data$reported_cases
plot(time, C, type = "l")


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



low2 <- loess(Y ~ X)
summary(low2)


############################

# polynomial
reg1 <- lm(Y ~ X)
rsquared <- summary(reg1)$r.squared 

reg2 <- lm(Y ~ X+ I(X^2))
rsquared <- c(rsquared, summary(reg2)$r.squared)

reg3 <- lm(Y ~ X + I(X^2) + I(X^3)) 
rsquared <- c(rsquared, summary(reg3)$r.squared)

reg4 <- lm(Y ~ X + I(X^2) + I(X^3) + I(X^4)) 
rsquared <- c(rsquared, summary(reg4)$r.squared)

reg5 <- lm(Y ~ X + I(X^2) + I(X^3) + I(X^4) + I(X^5)) 
rsquared <- c(rsquared, summary(reg5)$r.squared)

reg6 <- lm(Y ~ X + I(X^2) + I(X^3) + I(X^4) + I(X^5) + I(X^6)) 
rsquared <- c(rsquared, summary(reg6)$r.squared)

reg7 <- lm(Y ~ X + I(X^2) + I(X^3) + I(X^4) + I(X^5) + I(X^6) + I(X^7)) 
rsquared <- c(rsquared, summary(reg7)$r.squared)

reg8 <- lm(Y ~ X + I(X^2) + I(X^3) + I(X^4) + I(X^5) + I(X^6) + I(X^7) + I(X^8)) 
rsquared <- c(rsquared, summary(reg8)$r.squared)

reg9 <- lm(Y ~ X + I(X^2) + I(X^3) + I(X^4) + I(X^5) + I(X^6) + I(X^7) + I(X^8) + I(X^9)) 
rsquared <- c(rsquared, summary(reg9)$r.squared)

reg10 <- lm(Y ~ X + I(X^2) + I(X^3) + I(X^4) + I(X^5) + I(X^6) + I(X^7) + I(X^8) + I(X^9) + I(X^10)) 
rsquared <- c(rsquared, summary(reg10)$r.squared)

reg <- as.matrix(cbind(predict(reg1), predict(reg2), predict(reg3), predict(reg4), 
                       predict(reg5), predict(reg6), predict(reg7), predict(reg8), predict(reg9), predict(reg10)))
         
rho.all <- as.matrix(cbind(derivative(reg1, X), derivative(reg2, X), derivative(reg3, X), derivative(reg4, X),
                       derivative(reg5, X), derivative(reg6, X), derivative(reg7, X), derivative(reg8, X),
                       derivative(reg9, X), derivative(reg10, X)))       



score <- rep(NA, 10)
#penalty <- 5 * 10^-4 ### Check with Quentin
penalty <- 5*10^(-4)
for (n in 1:10) {
  score[n] <- rsquared[n] - (penalty*n^2)
}


maxn <- which.max(score)
plot(score, type = "l")
abline(v = maxn)



Yhat <- reg[ ,maxn]
rho <- rho.all[ ,maxn]
plot(1/rho, type = "l")


Z <- -(Yhat - Y)
plot(time, Z , type = "l", xlab = "Time", ylab = "Residual") #plot Z

# Ic are actual cases: reported cases multiplied by rho
Ic <- C*rho 
plot(Ic, type = "l")
lines(C, type = "l", col = "red")

plot(time, Y, type = "l", ylab = "Cumulative Births/Cases")
lines(time, Yhat, type = "l", col = "blue") #predicted cumulative births (blue)
lines(time, cumsum(C), col = "red")


## Keep only positive values
keep.index <- which(Ic[1:lengthdata - 1] > 0 & Ic[2:lengthdata] > 0)
length.pos <- length(keep.index)

## Set up the parameters
seas = rep(1:periodicity, lengthdata/periodicity)[1:lengthdata - 1][keep.index]
lInew <- log(Ic[2:lengthdata]) [keep.index]
lIold <- log(Ic[1:lengthdata -1]) [keep.index]
Zold = Z[1:lengthdata - 1] [keep.index]

sort(unique(seas))
plot(seas, type = "l")

######### This is from Bjornstadt (2007): Susceptible Reconstruction #######################
#Smean = seq(-0.015, 0.2, by=0.001)*.2E5  
min(Zold)
Smean = seq(abs(min(Zold)) + 1, abs(min(Zold))*3, by = 100) ## Reasonable range of candidate values
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
plot(Smean, -llik, type = "l", log = "y")
abline(v=Sbar)


lSold = log(Sbar + Zold) ## This is the best estimate (maximizes log likelihood)
# best estimates
lmfit <- lm(lInew ~ -1 + as.factor(seas) + lIold + offset(lSold))
sum <- summary(lmfit)
alpha <- lmfit$coef[periodicity + 1]

plot(exp(lIold), exp(lSold), cex = 0.6)

######################################################

#### Check other method: Taylor expansion
lmfit2 <- lm(lInew ~ -1 + as.factor(seas) + lIold + Zold)
sum2 <- summary(lmfit2)
plot(1:periodicity, exp(lmfit2$coef[1:periodicity]), type = "l", col = "red", 
     ylab = "Seasonal Coeff", xlab = "Week Number", main = "Iceland Seasonality")
lines(1:periodicity, exp(lmfit$coef[1:periodicity]) * Sbar, type = "l", lty = 2, col = "blue")
legend("bottomleft", legend = c("Using Sbar", "Using Z"), col = c("blue", "red"), lty = c(2,1), cex = 0.5)

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
plot(1:periodicity, lmfit$coef[1:periodicity] - mean(lmfit$coef[1:periodicity]), 
     type = "l", ylab = "seasonality", ylim = c(min(lower), max(upper)), 
     xaxp = c(0, periodicity+2, (periodicity+2)/2))
abline(h = 0)
lines(1:periodicity, upper, type = "l", lty = 2)
lines(1:periodicity, lower, type = "l", lty = 2)


##### Plotting the predicted I (from fitted model) vs actual I (corrected for underreporting)
keep2 <- Ic[1:lengthdata - 1] > 0 & Ic[2:lengthdata] > 0
keep <- which(Ic[2:lengthdata] > 0 & keep2 == TRUE)
lInew.pred <- rep(0, length(2:lengthdata))
lInew.pred[keep] <- predict(lmfit)

plot(time[1:lengthdata-1], Ic[1:lengthdata-1], type = "l", col = "red") #Observed cases corrected for underreporting
lines(time[2:lengthdata], exp(lInew.pred), type = "l") # Predicted cases


plot(exp(predict(lmfit)), exp(lInew), cex = 0.5, xlab = "Predicted Incidence", ylab = "Actual Incidence")
abline(a=0, b=1, col = "red")

################## Simulation #######################################

### Identify epidemics
# identify the starting index for each epidemic
if (sum(C <= sensitivity) / lengthdata > 0.2 ) {
  start.index <- ifelse(keep.index[2:length.pos] - keep.index[1:length.pos - 1] > 1, keep.index[2:length.pos], NA)
  start.index[1] <- keep.index[1]
  start.index <- start.index[is.na(start.index) == FALSE]
  end.index <- ifelse(keep.index[2:length.pos] - keep.index[1:length.pos - 1] > 1, keep.index[1:length.pos - 1], NA)
  end.index <- end.index[is.na(end.index) == FALSE]
  end.index <- append(end.index, keep.index[length(keep.index)])
  duration <- end.index - start.index
}

#create a matrix for each epidemic with start index, end index and duration
epi <- as.data.frame(as.matrix(cbind(start.index, end.index, duration)))

#Only keep epidemics which have a duration greater than 1 or more????
#epi <- epi[epi$duration > 1, ]


plot(Ic, type = "l")
for (i in 1:length(epi$start.index)) {
  abline(v = epi$start.index[i], col = "red")
}



# logI and logS are the full time series vectors
logI <- rep(0, lengthdata)
logI[keep.index] <- log(Ic[keep.index])
logS <- rep(0, lengthdata)
logS[keep.index] = log(Sbar + Z[keep.index])

# Set up parameters
predS <- matrix(NA, lengthdata, numsim)
predI <- matrix(NA, lengthdata, numsim)
lambda <- matrix(NA, lengthdata, numsim)
rstar <- exp(lrstar)



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


I <- rowMeans(predI)
S <- rowMeans(predS)
I <- ifelse(is.na(I) == TRUE, 0, I)
plot(I, type = "l", col = "red")
lines(Ic, type = "l", col = "blue")


plot(predI[,1], type = "l", col = "red", ylim = c(0, 3500))
for (i in 2:ncol(predI)) {
  lines(predI[,i], type = "l", col = "red")
}

lines(Ic, type = "l", col = "blue")



### Real vs predicted duration
pred.dur <- matrix(NA, length(epi$start.index), numsim)
pred.endindex <- matrix(NA, length(epi$start.index), numsim)

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
  }
  
  
}


pred.dur <- as.data.frame(pred.dur)
pred.dur$real.dur <- epi$duration
pred.dur$id <- c(1:length(epi$start.index))

duration <- reshape(pred.dur, 
               varying = list(names(pred.dur)[1:numsim]), 
               v.names = "pred.dur",
               idvar = c("id", "real.dur") ,     
               direction = "long")


# Real duration versus predicted duration
fit.dur <- lm(duration$pred.dur ~ -1 + duration$real.dur)
summary(fit.dur)
slope.dur <- round(fit.dur$coeff[1], 5)

plot(duration$real.dur, duration$pred.dur, pch = 1, cex = 0.5, 
     xlab = "Real Duration", ylab = "Predicted Duration", 
     main = paste("Slope = ", slope.dur, sep = ""), cex.main = 0.9)
abline(fit.dur, col = "blue")



### Real vs predicted size of each epidemic
# real size
for (e in 1:length(epi$start.index)) {
  epi$size[e] <- sum(Ic[epi$start.index[e]: epi$end.index[e]])
}


# predicted size
pred.size <- matrix(NA, length(epi$start.index), numsim)

for (c in 1:ncol(predI)) {
  predIcol <- predI[,c]
  for (s in 1:length(start.index.sim)){
    pred.size[s,c] <- sum(predIcol[epi$start.index[s]: pred.endindex[s,c]])
    }
}
  

pred.size <- as.data.frame(pred.size)
pred.size$real.size <- epi$size
pred.size$id <- c(1:length(epi$start.index))


epi$S0 <- exp(logS[epi$start.index])
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

plot(size$real.size, size$pred.size, pch = 1, cex = 0.5,
     xlab = "Real Size", ylab = "Predicted Size", 
     main = paste("Slope = ", slope.size, sep = ""), cex.main = 0.9)
abline(fit.size, col = "blue")

# S0 versus real size
fitS0 <- lm(epi$S0 ~ epi$size)
slope.S0.realsize <- round(fitS0$coeff[2], 5)
plot(epi$size, epi$S0, pch = 1, cex = 0.5,
     xlab = "Real Size", ylab = "S0", 
     main = paste("Slope = ", slope.S0.realsize, sep = ""), cex.main = 0.9)
abline(fitS0, col = "blue")

#S0 versus predicted size
fitS0.pred <- lm(size$realS0 ~ size$pred.size)
slope.S0.predsize <- round(fitS0.pred$coeff[2], 5)
plot(size$pred.size, size$realS0, pch = 1, cex = 0.5,
     xlab = "Predicted Size", ylab = "S0", 
     main = paste("Slope = ", slope.S0.predsize, sep = ""), cex.main = 0.9)
abline(fitS0.pred, col = "blue")


