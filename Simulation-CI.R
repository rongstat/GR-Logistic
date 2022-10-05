library(glmnet)
library(MASS)
f = function(x){
  exp(x)/(1+exp(x))
}

fp = function(x){
  exp(x)/(1+exp(x))^2
}

#sp=20, 25, 30, 35
sp=25

n=300

p.set = c(700,800,900,1000)
for(p.choose in 1:4){

  p=p.set[p.choose]


sig1 = toeplitz(seq(0.3, 0,length.out = p/10))
Sig = bdiag(rep(list(sig1),10))+diag(rep(0.7,p))


coverage.B=c()
coverage.Q=c()
coverage.R=c()
len.B=c()
len.Q=c()
len.R=c()
c.B=c()
c.Q=c()
c.R=c()
time.prop=c()

for(round in 1:100){

  b1 = rep(0,p)
  b1[1:sp] = runif(sp,-1,1)
  b2 = rep(0,p)
  b2[1:sp] = runif(sp,-1,1)
  
  
  
  (TrueB= as.numeric(b1 %*% Sig %*% b2))
  (TrueQ1= as.numeric(b1 %*% Sig %*% b1))
  (TrueQ2= as.numeric(b2 %*% Sig %*% b2))
  (TrueR= TrueB/sqrt(TrueQ1*TrueQ2))



X1=mvrnorm(n, mu=rep(0,p), Sigma=Sig)
prob = exp(X1 %*% b1)/(1+exp(X1 %*% b1))
y1 = rep(1,n)
while(sum(y1)/n<0.02 | sum(y1)/n>0.98 ){
  for(gen.y in 1:n){
    y1[gen.y]=rbinom(1,1,prob[gen.y])
  }
}



X2=mvrnorm(n, mu=rep(0,p), Sigma=Sig)
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


#################

coverage.B[round] = (B.hat-rho.B<TrueB) & (TrueB<B.hat+rho.B)
c.B[round] = B.hat
len.B[round] = 2*rho.B

coverage.Q[round] = (Q1.hat-rho.Q1<TrueQ1) & (TrueQ1<Q1.hat+rho.Q1)
c.Q[round] = Q1.hat
len.Q[round] = 2*rho.Q1

(coverage.R[round] = (R.hat-rho.R<TrueR) & (TrueR<R.hat+rho.R))
c.R[round] = R.hat
len.R[round] = min(c(1,R.hat+rho.R))-max(c(R.hat-rho.R,0))

}

print(round(c(mean(coverage.B),mean(len.B), mean(coverage.Q), mean(len.Q),mean(coverage.R),mean(len.R)),3))

}

