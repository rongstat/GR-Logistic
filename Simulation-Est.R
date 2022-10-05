args= commandArgs(trailingOnly=TRUE)


library(MASS)
library(CVXR)
library(AER)
library(Matrix);
library(glmnet);
library(flare);
library(hdi)

f = function(x){
  exp(x)/(1+exp(x))
}

sp=25
n=200
p=700



#sig1 = toeplitz(seq(0.3, 0,length.out = p/10))
#Sig = bdiag(rep(list(sig1),10))+diag(rep(0.7,p))

Sig = matrix(0.2,ncol=p,nrow=p)
diag(Sig)=1

X1=mvrnorm(n, mu=rep(0,p), Sigma=Sig)
X2=mvrnorm(n, mu=rep(0,p), Sigma=Sig)


out=matrix(nrow=500,ncol=12)
for(i in 1:500){
  #Set the regression coefficients.
  b1 = rep(0,p)
  b1[1:sp] = runif(sp,-1,1)
  b2 = rep(0,p)
  b2[1:sp] = runif(sp,-1,1)
  
  (TrueB= as.numeric(b1 %*% Sig %*% b2))
  (TrueQ1= as.numeric(b1 %*% Sig %*% b1))
  (TrueQ2= as.numeric(b2 %*% Sig %*% b2))
  (TrueR= TrueB/sqrt(TrueQ1*TrueQ2))
  
  
  
  prob = exp(X1 %*% b1)/(1+exp(X1 %*% b1))
  y1 = rep(1,n)
  while(sum(y1)/n<0.02 | sum(y1)/n>0.98 ){
    for(gen.y in 1:n){
      y1[gen.y]=rbinom(1,1,prob[gen.y])
    }
  }
  
prob = exp(X2 %*% b2)/(1+exp(X2 %*% b2))
y2 = rep(1,n)
while(sum(y2)/n<0.02 | sum(y2)/n>0.98 ){
  for(gen.y in 1:n){
    y2[gen.y]=rbinom(1,1,prob[gen.y])
  }
}



####### lasso estimation
fit = glmnet(x = X1, y = y1, family = "binomial", alpha = 1,  lambda = 0.12*sqrt(log(p)/n), intercept=F, standardize=F)
b.hat1 = as.vector(coef(fit))
fit = glmnet(x = X2, y = y2, family = "binomial", alpha = 1,  lambda = 0.12*sqrt(log(p)/n), intercept=F, standardize=F)
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


################### hdi1: be careful! This part should take Long time. Better run on separate cores.


fit.lasso = lasso.proj(x = X1, y = y1, standardize = F, family = "binomial")
b.hat1 = fit.lasso$bhat
fit.lasso = lasso.proj(x = X2, y = y2, standardize = F, family = "binomial")
b.hat2 = fit.lasso$bhat


B.db1 = b.hat1 %*% Sig.hat %*% b.hat2
Q1.db1 =  b.hat1 %*% Sig.hat %*% b.hat1
Q2.db1 =  b.hat2 %*% Sig.hat %*% b.hat2
R.db1 = B.db1/sqrt(Q1.db1*Q2.db1)
if(abs(R.db1)>1) R.db1=sign(R.db1)
if(Q1.db1*Q2.db1==0) R.db1=0


#####################hdi2: be careful! This part should take Long time. Better run on separate cores.
b.hat1 = ridge.proj(X1,y1, family="binomial", standardize = F)
b.hat2 = ridge.proj(X2,y2, family="binomial", standardize = F)
B.db2 = b.hat1$betahat %*% Sig.hat %*% b.hat2$betahat
Q1.db2 =  b.hat1$betahat %*% Sig.hat %*% b.hat1$betahat
Q2.db2 =  b.hat2$betahat %*% Sig.hat %*% b.hat2$betahat
R.db2 = B.db2/sqrt(Q1.db2*Q2.db2)
if(Q1.db2*Q2.db2==0){R.db2=0}else{
if(abs(R.db2)>1) R.db2=sign(R.db2)
}

############ plg

fit = cv.glmnet(x = X1, y = y1, family = "binomial", alpha = 1,  intercept=F, standardize=F)
b.hat1 = as.vector(coef(fit, s = "lambda.min"))
fit = cv.glmnet(x = X2, y = y2, family = "binomial", alpha = 1,  intercept=F, standardize=F)
b.hat2 = as.vector(coef(fit, s = "lambda.min"))
B.temp = t(as.matrix(b.hat1))[-1] %*% Sig.hat %*% (as.matrix(b.hat2))[-1]
Q1.temp = t(as.matrix(b.hat1))[-1] %*% Sig.hat %*% (as.matrix(b.hat1))[-1]
Q2.temp = t(as.matrix(b.hat2))[-1] %*% Sig.hat %*% (as.matrix(b.hat2))[-1]
############
Err.B = c(B.hat-TrueB, B.temp-TrueB, B.db1-TrueB, B.db2-TrueB)
Err.Q = c(Q1.hat-TrueQ1, Q1.temp-TrueQ1, Q1.1-TrueQ1, Q1.db2-TrueQ1)
Err.R = c(R.hat-TrueR, B.temp/sqrt(Q1.temp*Q2.temp)-TrueR, R.db1-TrueR, R.db2-TrueR)


(out[i,] = c(Err.B, Err.Q, Err.R))
}

colMeans(abs(out), na.rm=T)

save(out, file= args[1])
