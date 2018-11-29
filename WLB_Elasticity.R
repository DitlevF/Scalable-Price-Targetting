
rm(list = ls())


setwd('XX')

df <- readRDS(file = 'Seminar_BE/data.Rda')

library('gamlr')
library('glmnet')
library('dplyr')

# Estimation
colnames(df)[4:53] <- paste0('X',1:50)
X_cols <- colnames(df)[4:53]
price_vec <- c(19,39,59,79,99,159,199,249,299,399)
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
colnames(X)
x <- X[,X_cols] # For the elasticity.
y <- df[,'y']
dU <- df[,'dU']
x_bar <- colMeans(x)

df <- as.data.frame(df)

B <- 100 # Number of bootstraps
N <- length(y)

alpha_theta <- matrix(data = NA , nrow = ncol(x), ncol = B)
beta_theta <- matrix(data = NA , nrow = ncol(x), ncol = B)
lambda <- rep(NA,B)

elasticity_discrete <- matrix(data = NA, nrow = B, ncol = length(price_vec))
elasticity_cont <- matrix(data = NA, nrow = B, ncol = length(price_vec))
Average_elasticity <- matrix(data = NA, nrow = B, ncol = 1)


for(i in 1:B){
  cat("iter =", i, "\r")
  
  V <- rexp(N, rate = 1) 
  df_bs <- sample_n(df, N, replace = TRUE, weight = V) #Draw Weights
  X_bs <- subset(df_bs, select = -c(y, prices,dU)) 
  y_bs <- df_bs[,'y']
  
  fit <- cv.gamlr(X_bs,y_bs, family = "binomial", intercept = FALSE, verb = FALSE) #run gamlr
  
  alpha_theta[,i] <- coef(fit)[2:51]
  beta_theta[,i] <-  coef(fit)[52:101]
  lambda[i] <- fit$lambda.1se
  
  for(j in 1:length(price_vec)){
    
    dU_0 <- x %*% alpha_theta[,i] + (x %*% beta_theta[,i]) * price_vec[j]
    dU_1 <- x %*% alpha_theta[,i] + (x %*% beta_theta[,i]) * (price_vec[j]*1.01)
    
    elasticity_discrete[i,j] <- (((sum(as.numeric(dU_1>0)) - sum(as.numeric(dU_0)>0))) / sum(as.numeric(dU_0)>0)) * 100
  }
}

theta <- rbind(alpha_theta,beta_theta)

saveRDS(theta, file = 'thetas.Rda')
saveRDS(elasticity_discrete, file = 'elasticities.Rda')
saveRDS(lambda, file = 'lambdas.Rda')


