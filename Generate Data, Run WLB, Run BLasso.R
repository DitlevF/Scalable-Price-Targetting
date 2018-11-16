# Bayesian Lasso Elasticity
rm(list = ls())

setwd("Set your WD here")
set.seed(251093)

library('coda')
library('gamlr')
library('glmnet')
library('dplyr')
library('tmvtnorm')


###Simulate Data###
n <- 10000
k <- 50
sigma_zeta <- 60

# Willingness-to-Pay (WTP) heterogeneity as a function of observable characteristics.

active_alpha <- c(1,0.5,0.5,0.4,0.7,0.2,0.2) # Better. 
alpha <- c(active_alpha, rep(0,5), active_alpha, rep(0,k-5-2*length(active_alpha)))


active_beta <- c(-0.002,-0.007,-0.005,-0.003) 
beta <- c(rep(0,5), active_beta,rep(0,k-5-2*length(active_beta)), active_beta) 

# Generate correlation matrix
corr_x <- matrix(nrow = k, ncol = k)
for(i in 1:k) for(j in 1:k) corr_x[i,j] <- 0.5**abs(i-j)


x_low <- rtmvnorm(0.65*n, mean = c(rep(1,5), rep(100,length(active_beta)), rep(1,k-5-2*length(active_beta)),
                                  rep(100,length(active_beta))), sigma = corr_x, lower = rep(0,k))
x1 <-rtmvnorm(0.1*n, mean = rep(1,k), sigma = corr_x, lower = rep(0,k))
x2 <-rtmvnorm(0.05*n, mean = rep(20,k), sigma = corr_x, lower = rep(0,k))
x3 <-rtmvnorm(0.2*n, mean = rep(50,k), sigma = corr_x, lower = rep(0,k))
x <- rbind(x_low,x1,x2,x3)


zeta <- rnorm(n, mean = 0, sd = sigma_zeta)


# Generate random price data for each observation
price_vec <- c(19,39,59,79,99,159,199,249,299,399)
blocks <- (1/length(price_vec))*1:length(price_vec)
u <- runif(n)
prices <- rep(NA,n)

for(i in 1:n){
  for(j in 1:length(price_vec)){
    if(u[i] - blocks[j] <= 0 ){
      prices[i] <- price_vec[j]
      break
    }
  }
}

dU <- x %*% alpha + (x %*% beta) * prices + zeta

conversion_rate <- matrix(data = NA, nrow = 1, ncol = length(price_vec ))
colnames(conversion_rate) <- c(as.character(price_vec))

for(i in 1:length(price_vec)){
  conversion_rate[1,i] <- (length(which(prices[which(dU > 0)] == price_vec[i])) / 
                                    length(which(prices == price_vec[i])))
 
}
conversion_rate

# Signal-to-noise ratio
StN <- sqrt(sum((x %*% alpha + x %*% beta) ^ 2)) / (sqrt(n) * sigma_zeta)

# Final (observed) data set
y <- as.integer(dU > 0)
data <- cbind(dU,y,prices)
colnames(data)[1] <- 'dU'
df <- as.data.frame(cbind(data,x))

# Estimation
colnames(df)[4:53] <- paste0('X',1:50)
X_cols <- colnames(df)[4:53]
price_vec <- 80:399

# Create price-covariate interaction
for(i in X_cols){
  new_col <- paste0(i,"Xprice")
  df <- df %>%
    mutate(z = (get(i))*prices) 
  
  names(df)[names(df) == 'z'] <- new_col
}


df <- as.matrix.data.frame(df)
X <- subset(df, select = -c(y, prices,dU))
x <- X[,X_cols] # For the elasticity.
X <- scale(X, scale = FALSE)
y <- df[,'y']
dU <- df[,'dU']
df <- as.data.frame(df)

N <- length(y)

### Bayesian Lasso ###
draws <- blasso_probit(X,y, iter = 30000, fix_sd = 1)

betas <- draws$betas

betas <- betas[28001:30000,]

post_betas <- as.mcmc(betas)
plot(post_betas[,50:70])
plot(post_betas[,c(5,24,51,85)])
plot(post_betas[,90:100])


# Calculating the elasticity
alpha_theta_lasso <- t(betas[,1:50])
beta_theta_lasso <- t(betas[,51:100])
no_draws <- dim(betas)[1]

elasticity_lasso <- matrix(data = NA, nrow = no_draws, ncol = length(price_vec))
profit_lasso <- matrix(data = NA, nrow = no_draws, ncol = length(price_vec))
price_opt_lasso <- rep(NA, no_draws)

for(i in 1:no_draws){
  cat("iter =", i, "\r")
    for(j in 1:length(price_vec)){
    
      dU_0 <- x %*% alpha_theta_lasso[,i] + (x %*% beta_theta_lasso[,i]) * price_vec[j]
      dU_1 <- x %*% alpha_theta_lasso[,i] + (x %*% beta_theta_lasso[,i]) * (price_vec[j]*1.01)
    
      elasticity_lasso[i,j] <- (((sum(as.numeric(dU_1>0)) - sum(as.numeric(dU_0)>0))) / sum(as.numeric(dU_0)>0)) * 100
      profit_lasso[i,j] <- sum((dU_0>0))*price_vec[j]
    }
  price_opt_lasso[i] <- price_vec[which.max(profit_lasso[i,])]
}


### Weighted Likelihood Bootstrap###

B <- 100 #Number of Bootstraps

#Generate Containers
alpha_theta <- matrix(data = NA , nrow = ncol(x), ncol = B)
beta_theta <- matrix(data = NA , nrow = ncol(x), ncol = B)
elasticity_discrete <- matrix(data = NA, nrow = B, ncol = length(price_vec))
Profit <- matrix(data = NA, nrow = B, ncol = length(price_vec))
Price_Opt <- rep(NA, B)

#Run WLB. Calculate elasticity, profit and profit maximising price.
for(i in 1:B){
  cat("iter =", i, "\r")
  V <- rexp(N, rate = 1) 
  df_bs <- sample_n(df, N, replace = TRUE, weight = V) #Draw Weights
  X_bs <- subset(df_bs, select = -c(y, prices,dU)) 
  y_bs <- df_bs[,'y']
  
  fit <- cv.gamlr(X_bs,y_bs, family = "binomial", intercept = FALSE, verb = FALSE) #run gamlr
  
  alpha_theta[,i] <- coef(fit)[2:51]
  beta_theta[,i] <-  coef(fit)[52:101]
  
  for(j in 1:length(price_vec)){
    #Calculate Profit and Elasticity at different prices
    dU_0 <- x %*% alpha_theta[,i] + (x %*% beta_theta[,i]) * price_vec[j]
    dU_1 <- x %*% alpha_theta[,i] + (x %*% beta_theta[,i]) * (price_vec[j]*1.01)
    
    elasticity_discrete[i,j] <- (sum(dU_1>0)-sum(dU_0>0)) / ((price_vec[j]*1.01) - price_vec[j]) * price_vec[j] / sum(dU_0>0)
    Profit[i,j] <- sum((dU_0>0))*price_vec[j]
  }
  Price_Opt[i] <- price_vec[which.max(Profit[i,])]
} 

###Calculate Elasticity using true data###

#Profit and elasticity given "real" alpha and beta
profit_real <- rep(NA, length(price_vec))
elasticity_real <- rep(NA, length(price_vec))
for(j in 1:length(price_vec)){
  dU_real <- x %*% alpha + (x %*% beta) * price_vec[j]
  dU_real2 <- x %*% alpha + (x %*% beta) * (price_vec[j]*1.01)
  elasticity_real[j] <- (sum(dU_real2>0)-sum(dU_real>0)) / ((price_vec[j]*1.01) - price_vec[j]) * price_vec[j] / sum(dU_real>0)
  profit_real[j] <- sum((dU_real>0))*price_vec[j]
}


###MSE, ELasticities at mean and profit at mean####

##Mean Squared Error##

#Generate Containers
SE_WLB     <- matrix(data=NA, nrow=B, ncol=length(price_vec))
MSE_WLB    <- rep(NA, length(price_vec))
SE_BLasso  <- matrix(data=NA, nrow=no_draws, ncol=length(price_vec)) 
MSE_BLasso <- rep(NA, length(price_vec))
for(j in 1:length(price_vec)){
  for(i in 1:B){
    SE_WLB[i,] <- (elasticity_real - elasticity_discrete[i,])**2
  }
  MSE_WLB[j] <- mean(SE_WLB[,j])
  for(k in 1:no_draws){
    SE_BLasso[k,] <- (elasticity_real - elasticity_lasso[k,])**2
  }
  MSE_BLasso[j] <- mean(SE_BLasso[,j])
}

#Total MSE
Total_MSE_WLB <- mean(MSE_WLB)
Total_MSE_WLB
Total_MSE_BLasso <- mean(MSE_BLasso)
Total_MSE_BLasso

#MSE at intervals
Interval_MSE_BLasso <- NA
Interval_MSE_WLB <- NA
for(i in c(50,100,150,200,250,300)){
  Interval_MSE_BLasso[i] <- mean(MSE_BLasso[(i-50):i])
  Interval_MSE_WLB[i]    <- mean(MSE_WLB[(i-50):i])
}
Interval_MSE_BLasso <- Interval_MSE_BLasso[!is.na(Interval_MSE_BLasso)]
Interval_MSE_BLasso
Interval_MSE_WLB <- Interval_MSE_WLB[!is.na(Interval_MSE_WLB)]
Interval_MSE_WLB


#Average of profit over WLB bootstraps
prof <- rep(NA,length(price_vec))
for(i in 1:length(price_vec)){
  prof[i] <- sum(as.numeric(Profit[,i]))/B
}

#Average of profit over BL draws
prof2 <- rep(NA,length(price_vec))
for(i in 1:length(price_vec)){
  prof2[i] <- sum(as.numeric(profit_lasso[,i]))/no_draws
}

#Average of elasiticies over WLB bootstraps
elas_WLB <- rep(NA,length(price_vec))
for(i in 1:length(price_vec)){
  elas_WLB[i] <- sum(elasticity_discrete[,i])/B
}

#Average of elasticities over BL Draws
elas_BLasso <- rep(NA,length(price_vec))
for(i in 1:length(price_vec)){
  elas_BLasso[i] <- sum(elasticity_lasso[,i])/no_draws
}

save.image("data.RData")





