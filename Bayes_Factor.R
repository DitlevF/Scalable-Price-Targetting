# Bayes Factors

rm(list = ls())
# Absolutely Continuous Spikes
library(coda)
library(dplyr)

# Load working directory either at home or at CSS
setwd('XX')

# Load Signal-to-Noise < 1
#load("Seminar_BE/Resultater/Bayes Factor Data/StN 10 - 10000 Obs.RData")

# Load Signal-to-Noise = 3
#load("Seminar_BE/Resultater/Forskellige Signal-to-noise/Ca. 3/SNR3.RData")

# Load Signal-to-Noise = 10
#load("Seminar_BE/Resultater/Forskellige Signal-to-noise/10/10.RData")


df <- as.matrix.data.frame(df)
X <- subset(df, select = -c(y, prices,dU)); X <- as.matrix(X)
y <- df[,'y']

logL_probit <- function(X,y, theta){
  XB <- X %*% theta
  FXB <- pnorm(XB,0,1)
  logL <- rep(NA,length(y))
  for(i in 1:length(y)){
    if(y[i] == 1){
      logL[i] <- log(FXB[i])
    } else{
      logL[i] <- log(1-FXB[i])
    }
  }
  L <- sum(logL)
  return(L)
}
logL_logit <- function(X,y, theta){
  XB <- X %*% theta
  FXB <- 1/(1+exp(-XB))
  logL <- rep(NA,length(y))
  for(i in 1:length(y)){
    if(y[i] == 1){
      logL[i] <- log(FXB[i])
    } else{
      logL[i] <- log(1-FXB[i])
    }
  }
  L <- sum(logL)
  return(L)
}

# WLB

thetas <- rbind(alpha_theta,beta_theta)
thetas <- t(as.matrix(thetas))
m_WLB <- dim(thetas)[1]

# Calculate log-likelihoods
logl_wlb <- rep(NA,m_WLB)
for(i in 1:m_WLB) logl_wlb[i] <- logL_logit(X,y,thetas[i,])
alpha <- min(logl_wlb)
l_logit_wlb <- exp(-logl_wlb + alpha)

### Bayesian LASSO ###
betas <- rbind(alpha_theta_lasso,beta_theta_lasso)
betas <- t(as.matrix(betas)); bs <- dim(betas)[1]
betas <- betas[(bs-1999):bs,]
m_BL <- dim(betas)[1]

# Calculate log-likelihoods
logl_bl <- rep(NA,m_BL)
for(i in 1:m_BL) logl_bl[i] <- logL_probit(X,y,betas[i,])
l_probit_blasso <- exp(-logl_bl + alpha)

BF <- (m_BL/m_WLB) * (sum(l_logit_wlb) / sum(l_probit_blasso))
log(BF)