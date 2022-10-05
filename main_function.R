

##### Input:
# X1: covariate matrix for data 1, with rows being samples and columns being features.
# y1: response vector for data 1
# X2: covariate matrix for data 2, with rows being samples and columns being features.
# y2: response vector for data 2

##### Output:
# g.corr, g.cov, g.var1, g.var2 are point estimators for genetic correlation, covariance, and genetic variances for data 1 and 2, respectively. 
# CI.g.corr, CI.g.cov, CI.g.var1, CI.g.var2 are 95% confidence intervals for genetic correlation, covariance, and genetic variances for data 1 and 2, respectively. 

genetic_correlation <- function(X1=X1, y1=y1, X2=X2, y2=y2, lambda.lasso = 0.12){
  
  ####### lasso estimation
  
  fit = glmnet(x = X1, y = y1, family = "binomial", alpha = 1,  lambda = lambda.lasso*sqrt(log(p)/n), intercept=F, standardize=F)
  b.hat1 = as.vector(coef(fit))
  fit = glmnet(x = X2, y = y2, family = "binomial", alpha = 1,  lambda = lambda.lasso*sqrt(log(p)/n), intercept=F, standardize=F)
  b.hat2 = as.vector(coef(fit))
  
  #Covariance estimation
  Sig.hat = cov( rbind(X1,X2))
  
  #Debiased estimators
  B.temp = t(as.matrix(b.hat1))[-1] %*% Sig.hat %*% (as.matrix(b.hat2))[-1]
  D1=mean(fp(X1%*% t(as.matrix(b.hat1))[-1]+b.hat1[1])^(-1) *(f(X1%*% t(as.matrix(b.hat1))[-1]+b.hat1[1])-y1)* X1 %*% (as.matrix(b.hat2))[-1]) 
  D2=mean(fp(X2%*% t(as.matrix(b.hat2))[-1]+b.hat2[1])^(-1) *(f(X2%*% t(as.matrix(b.hat2))[-1]+b.hat2[1])-y2)* X2 %*% (as.matrix(b.hat1))[-1]) 
  B.hat = B.temp-D1-D2
  
  Q1.temp = t(as.matrix(b.hat1))[-1] %*% Sig.hat %*% (as.matrix(b.hat1))[-1]
  D3=mean(fp(X1%*% t(as.matrix(b.hat1))[-1]+b.hat1[1])^(-1) *(f(X1%*% t(as.matrix(b.hat1))[-1]+b.hat1[1])-y1)* X1 %*% (as.matrix(b.hat1))[-1]) 
  Q1.hat = Q1.temp - 2*D3
  
  Q2.temp = t(as.matrix(b.hat2))[-1] %*% Sig.hat %*% (as.matrix(b.hat2))[-1]
  D4=mean(fp(X2%*% t(as.matrix(b.hat2))[-1]+b.hat2[1])^(-1) *(f(X2%*% t(as.matrix(b.hat2))[-1]+b.hat2[1])-y2)* X2 %*% (as.matrix(b.hat2))[-1]) 
  Q2.hat = Q2.temp - 2*D4
  
  
  R.hat = B.hat/sqrt(Q1.hat*Q2.hat)
  
  if(abs(R.hat)>1) R.hat=sign(R.hat)
  if(Q1.hat*Q2.hat==0) R.hat=0
  
  v.B = 2 * mean(fp(X1%*% t(as.matrix(b.hat1))[-1]+b.hat1[1])^(-1) * (X1 %*% t(as.matrix(b.hat1))[-1]) * (X1 %*% t(as.matrix(b.hat1))[-1])) + 
    2 * mean(fp(X2%*% t(as.matrix(b.hat2))[-1]+b.hat2[1])^(-1) * (X2 %*% t(as.matrix(b.hat2))[-1]) * (X2 %*% t(as.matrix(b.hat2))[-1])) +
    (sum(((X1 %*% t(as.matrix(b.hat1))[-1]) * (X1 %*% t(as.matrix(b.hat2))[-1])-
            c(t(as.matrix(b.hat1))[-1] %*% Sig.hat %*% (as.matrix(b.hat2))[-1]))^2) +
       sum(((X2 %*% t(as.matrix(b.hat1))[-1]) * (X2 %*% t(as.matrix(b.hat2))[-1])-
              c(t(as.matrix(b.hat1))[-1] %*% Sig.hat %*% (as.matrix(b.hat2))[-1]))^2))/2/n
  
  rho.B = 1.96 * sqrt(v.B) /sqrt(2*n)
  
  v.Q1 = 8*mean(fp(X1%*% t(as.matrix(b.hat1))[-1]+b.hat1[1])^(-1) * (X1 %*% t(as.matrix(b.hat1))[-1]) * (X1 %*% t(as.matrix(b.hat1))[-1]))  + 
    (sum(((X1 %*% t(as.matrix(b.hat1))[-1]) * (X1 %*% t(as.matrix(b.hat1))[-1])-
            c(t(as.matrix(b.hat1))[-1] %*% Sig.hat %*% (as.matrix(b.hat1))[-1]))^2) +
       sum(((X2 %*% t(as.matrix(b.hat1))[-1]) * (X2 %*% t(as.matrix(b.hat1))[-1])-
              c(t(as.matrix(b.hat1))[-1] %*% Sig.hat %*% (as.matrix(b.hat1))[-1]))^2))/2/n
  
  
  rho.Q1 = 1.96 * sqrt(v.Q1) /sqrt(2*n)
  
  rho.R = rho.B/sqrt(Q1.hat*Q2.hat)
  
  return(list(g.corr = R.hat, g.cov = B.hat, g.var1 = Q1.hat, g.var2 = Q2.hat,
              CI.g.corr = c(R.hat-rho.R, R.hat+rho.R), 
              CI.g.cov = c(B.hat-rho.B, B.hat+rho.B),
              CI.g.var1 = c(Q1.hat-rho.Q1, Q1.hat+rho.Q1),
              CI.g.var2 = c(Q2.hat-rho.Q2, Q2.hat+rho.Q2))
         )
  
  
}