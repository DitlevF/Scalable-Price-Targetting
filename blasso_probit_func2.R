blasso_probit <- function(X,y, fix_sd = FALSE, r = 1, delta = 1.78, iter = 10000){
  require(truncnorm)
  require(VGAM)
  
  # Probit extension
  draw_ystar <- function(X, Y, alpha, mu, sd = 1){
    n <- length(Y);
    mean_XB <- X %*% alpha + mu
    K1 <- length(Y[Y == 1]); K2 <- length(Y[Y == 0])
    y_star <- rep(0,n);
    
    a1 <- rep(0,K1); b1 <- rep(Inf,K1); a2 <- rep(-Inf,K2); b2 <- rep(0,K2); sd1 <- rep(sd, K1); sd2 <- rep(sd, K2)
    y_star[Y == 1] <- rtruncnorm(1, a = a1, b = b1, mean = mean_XB[Y == 1], sd = sd1)
    y_star[Y == 0] <- rtruncnorm(1, a = a2, b = b2, mean = mean_XB[Y == 0], sd = sd2)
    return(y_star)
  }
  
  # Parameters
  N <- dim(X)[1]
  k <- dim(X)[2]
  XtX <- t(X) %*% X
  
  # Storage Matrices
  beta_draws <- matrix(data = NA, nrow = iter, ncol = k)
  sigma_draws <- rep(NA, iter)
  tau_draws <- matrix(data = NA, nrow = iter, ncol = k)
  lambda_draws <- rep(NA, iter)
  
  # Initialize
  lambda2 <- 1; mu_dot_init <- 1 
  tau2 <- 1/rinv.gaussian(k, mu_dot_init, lambda2)
  D <- diag(tau2)
  sigma2 <- 1
  pr_fit <- glm.fit(X,y,family = binomial(link = "probit")); beta <- pr_fit$coefficients # Initialize betas at proper values
  mu <- 0
  #beta <- runif(k,-0.5,0.5)
  
  
  for(i in 1:iter){
    
    if (i %% 100 == 0) {
      cat('Iteration:', i, "\r")
    }
    
    # Probit Step
    Y_star <- if(fix_sd == FALSE) draw_ystar(X,y,beta, mu, sd = 1) else draw_ystar(X,y,beta, mu, sd = fix_sd)
    yc <- Y_star - mean(Y_star)
    
    # Sample beta
    A <- XtX + solve(D)
    beta <- c(rmvnorm(1, mean = solve(A) %*% t(X) %*% yc, sigma = sigma2 * solve(A) ))
    beta_draws[i,] <- beta
    
    # Sample sigma2
    
    if(fix_sd == FALSE){
      sh <- (N-1)/2 + k/2
      sc <- 0.5*t(yc - X %*% beta) %*% (yc - X %*% beta) + 0.5*t(beta) %*% solve(D) %*% beta
      sigma2 <- 1/rgamma(1, sh, sc )
    } else{
      sigma2 <- fix_sd^2
    }
    
    sigma_draws[i] <- sigma2
    
    # Sample tau2
    mu_dot <- sqrt(lambda2 * sigma2 / beta^2)
    tau2 <- 1/rinv.gaussian(k, mu_dot, lambda2)
    
    # Sample lambda2
    lambda2 <- rgamma(1, shape = k + r, rate = 0.5 * sum(tau2) + delta)
    lambda_draws[i] <- lambda2

    # Sample mu
    mu <- rnorm(1,mean(Y_star), sigma2/N)
  }
  
  draws <- list('betas' = data.frame('betas' = beta_draws), 'sigma2' = data.frame('sigma2' = sigma_draws), 
                'lambda' = lambda_draws)
  
  colnames(draws$betas) <- colnames(X)
  return(draws)
}