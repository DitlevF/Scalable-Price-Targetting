# Generate artificial ZipRecruiter data
rm(list = ls())

setwd("Set your own WD here")


set.seed(251093)
library('tmvtnorm')
n <- 10000
k <- 50
sigma_zeta <- 60
#sigma_zeta <- pi/3

# Willingness-to-Pay (WTP) heterogeneity as a function of observable characteristics.

active_alpha <- c(250,250,250,320,200,100,100)
active_alpha <- c(1,0.5,0.5,0.4,0.7,0.2,0.2) # Better. 

alpha <- c(active_alpha, rep(0,5), active_alpha, rep(0,k-5-2*length(active_alpha)))

active_beta <- c(-2.0,-1.5,-1.5,-2)
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
#print(xtable(conversion_rate, type = "latex"), file = "filename2.tex")

# Signal-to-noise ratio
sqrt(sum((x %*% alpha + x %*% beta) ^ 2)) / (sqrt(n) * sigma_zeta)

# Final (observed) data set
y <- as.integer(dU > 0)
data <- cbind(dU,y,prices)
colnames(data)[1] <- 'dU'
data <- as.data.frame(cbind(data,x))

save.image("Set WD here")

