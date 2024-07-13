require(pscl)
library(glmnet)
library(MASS)
library(MCMCpack)

n <- 300  # the number of observations
p <- 5  # the number of predictor variables
om1 <- 0.1    # zero-inflation prob
om2 <- 0.05    # large outlier prob
aa <- 20    # outlier location (10 or 20)


# settings
Beta <- rep(0,p)  
Beta[c(1,3,5)] <- c(-0.2, 0.3, 0.6)   # non-zero coefficients
int <- 2   # intercept
Para <- c(int, Beta)


# covariates
rho <- 0.2
mat <- matrix(NA, p, p)
for(k in 1:p){
  for(j in 1:p){ mat[k,j]=rho^(abs(k-j)) }
}

X <- mvrnorm(n, rep(0,p), mat) 


# response values
Mu <- as.vector( exp(int + X%*%Beta) )
Y <- rpois(n, Mu)

plot(X[,3], Y)

# outlier and zero inflation
ch <- sample(c(0, 1, 2), size=n, prob=c(1-om1-om2, om1, om2), replace=T)
Y[ch==1] <- 0
Y[ch==2] <- Y[ch==2] + aa

plot(X[,3], Y)





##approach algorithm
lasso_possion_EM <- function(X, Y){
  ## EM 
  dEH <- function(x){
    (log(1+x))^(-1/2)/(1+x)/(1+log(1+x))/3.14
  }
  
  mc <- 1000
  
  
  
  cv_fit_I <- cv.glmnet(X, Y ,family="poisson", alpha=1)
  best_lambda_I <- cv_fit_I$lambda.min
  final_fit_I <- glmnet(X, Y,family="poisson", alpha=1, lambda = best_lambda_I)
  hBeta <- as.vector(coef(final_fit_I))
  
  
  ss <- 0.05
  maxit <- 100
  n <- length(Y)
  
  for(j in 1:maxit){
    hBeta.old <- hBeta
    # E-step
    lam <- exp(as.vector(cbind(1,X)%*%hBeta))
    mD1 <- c()
    mD2 <- c()
    for(i in 1:n){
      rn <- rgamma(mc, 1+Y[i], lam[i])
      mD1[i] <- mean(dEH(rn))/lam[i]
      mD2[i] <- mean(rn*dEH(rn))/lam[i]
    }
    mZ <- ss*mD1 / ((1-ss)*dpois(Y, lam) + ss*mD1)
    mEta <- mD2/mD1
    
    # M-step
    ss <- mean(mZ)
    offset <- log(1-mZ+mZ*mEta)
    hBeta_o <- coef(glm(Y~X, offset=offset, family="poisson"))
    
    #with lasso  
    cv_fit <- cv.glmnet(X, Y, offset=offset ,family="poisson", alpha=1)
    best_lambda <- cv_fit$lambda.min
    final_fit <- glmnet(X, Y, offset=offset ,family="poisson", alpha=1, lambda = best_lambda)
    
    hBeta <- as.vector(coef(final_fit))
    
    # convergence check 
    dd <- 100*sum(abs(hBeta-hBeta.old)) / sum(abs(hBeta)+0.0001)
    # print(hBeta)
    if(dd<0.005){ break }
  }
  outlier_idx <- ifelse(mZ >0.5, 1,0)
  
  
  return(list(beta = hBeta, beta_without_lasso = hBeta_o, Z_latent = mZ, outlier_idx = outlier_idx, best_M = final_fit, eta = mEta))
  
}






#calculate DM
lasso_possion_EM_Y <- function(X, Y){
  dEH <- function(x){
    (log(1+x))^(-1/2)/(1+x)/(1+log(1+x))/3.14
  }
  
  mc <- 1000
  hBeta <- coef(glm(Y~X, family="poisson"))
  ss <- 0.05
  maxit <- 1
  n <- length(Y)
  
  for(j in 1:maxit){
    hBeta.old <- hBeta
    # E-step
    lam <- exp(as.vector(cbind(1,X)%*%hBeta))
    mD1 <- c()
    mD2 <- c()
    for(i in 1:n){
      rn <- rgamma(mc, 1+Y[i], lam[i])
      mD1[i] <- mean(dEH(rn))/lam[i]
      mD2[i] <- mean(rn*dEH(rn))/lam[i]
    }
    mZ <- ss*mD1 / ((1-ss)*dpois(Y, lam) + ss*mD1)
    mEta <- mD2/mD1
    
    # M-step
    ss <- mean(mZ)
    offset <- log(1-mZ+mZ*mEta)
    hBeta_o <- coef(glm(Y~X, offset=offset, family="poisson"))
    
    #with lasso  
    cv_fit <- cv.glmnet(X, Y, offset=offset ,family="poisson", alpha=1)
    best_lambda <- cv_fit$lambda.min
    final_fit <- glmnet(X, Y, offset=offset ,family="poisson", alpha=1, lambda = best_lambda)
    
    hBeta <- as.vector(coef(final_fit))
    
    # convergence check 
    dd <- 100*sum(abs(hBeta-hBeta.old)) / sum(abs(hBeta)+0.0001)
    # print(hBeta)
    if(dd<0.005){ break }
  }
  outlier_idx <- ifelse(mZ >0.5, 1,0)
  
  
  return(list(beta = hBeta, beta_without_lasso = hBeta_o, Z_latent = mZ, outlier_idx = outlier_idx, best_M = final_fit, eta = mEta))
  
}



r_1 <- lasso_possion_EM(X,Y)


lasso_possion_EM_DM_Y <- function(X, Y, hBeta){
  dEH <- function(x){
    (log(1+x))^(-1/2)/(1+x)/(1+log(1+x))/3.14
  }
  
  mc <- 1000
  ss <- 0.05
  maxit <- 1
  n <- length(Y)
  
  for(j in 1:maxit){
    hBeta.old <- hBeta
    # E-step
    lam <- exp(as.vector(cbind(1,X)%*%hBeta))
    mD1 <- c()
    mD2 <- c()
    for(i in 1:n){
      rn <- rgamma(mc, 1+Y[i], lam[i])
      mD1[i] <- mean(dEH(rn))/lam[i]
      mD2[i] <- mean(rn*dEH(rn))/lam[i]
    }
    mZ <- ss*mD1 / ((1-ss)*dpois(Y, lam) + ss*mD1)
    mEta <- mD2/mD1
    
    # M-step
    ss <- mean(mZ)
    offset <- log(1-mZ+mZ*mEta)
    hBeta_o <- coef(glm(Y~X, offset=offset, family="poisson"))
    
    #with lasso  
    cv_fit <- cv.glmnet(X, Y, offset=offset ,family="poisson", alpha=1)
    best_lambda <- cv_fit$lambda.min
    final_fit <- glmnet(X, Y, offset=offset ,family="poisson", alpha=1, lambda = best_lambda)
    
    hBeta <- as.vector(coef(final_fit))
    
    # convergence check 
    dd <- 100*sum(abs(hBeta-hBeta.old)) / sum(abs(hBeta)+0.0001)
    # print(hBeta)
    if(dd<0.005){ break }
  }
  outlier_idx <- ifelse(mZ >0.5, 1,0)
  
  
  return(list(beta = hBeta, beta_without_lasso = hBeta_o, Z_latent = mZ, outlier_idx = outlier_idx, best_M = final_fit, eta = mEta))
  
}














I_beta <- function(beta, omega, eta, x, y) {
  n <- length(y)
  p <- length(beta)
  lambda <- exp(x %*% beta)
  weight <- omega * eta + (1 - omega)
  
  I_beta_mat <- matrix(0, p, p)
  for (i in 1:n) {
    I_beta_mat <- I_beta_mat + weight[i] * lambda[i] * x[i,] %*% t(x[i,])
  }
  
  return(I_beta_mat)
}





DM <- function(beta, omega, eta, x, y, eps=1e-6) {
  p <- length(beta)
  DM_mat <- matrix(0, p, p)
  
  
  
  for (j in 1:p) {
    beta_plus <- beta
    beta_plus[j] <- beta_plus[j] + eps
    theta_plus <- c(beta_plus, omega, eta)
    theta_new <- lasso_possion_EM_DM_Y(X = X, Y = y, beta_plus)
    DM_mat[j,] <- c(theta_new$beta_without_lasso - beta) / eps
  }
  
  return(DM_mat)
}


C_I <- function(mu , sd, n){
  z <- qnorm(0.975)  # 95% 置信度对应的 z 值
  lower_ci <- mu - z * sd / sqrt(n)
  upper_ci <- mu + z * sd / sqrt(n)
  return(c(lower_ci,upper_ci))
}

v_m <- function(beta,eta,Z_latent,X,Y){
  omega <- mean(Z_latent)
  f_x <- cbind(1,X)
  I_c_mat <- I_beta(beta = beta,omega = omega,eta = eta,x = f_x, y = Y)
  DM_mat <- DM(beta = beta ,omega = omega, eta = eta ,x = X, y = Y)
  
  I_c_inv <- solve(I_c_mat)
  matrix_v <- I_c_inv + I_c_inv %*% DM_mat %*% solve(diag(p + 1) - DM_mat) %*% I_c_inv
  
  s_d <- sqrt(diag(matrix_v))
  CI <-c()
  for(i in 1:length(beta)){
    CI <- cbind(CI, C_I(beta[i], s_d[i],300))
  }
  return(list(M_p = sqrt(diag(matrix_v)),ci = CI))
  
}



beta_f <- r_1$beta_without_lasso
v_m(beta = r_1$beta_without_lasso,eta = r_1$eta,Z_latent = r_1$Z_latent, X = X, Y = Y)
