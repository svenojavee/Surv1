rm(list=ls())
gc()

library(ars)
library(invgamma)
library(data.table)
library(MASS)

#Load the ARS functions
source("/Users/admin/Desktop/Surv1/ARS_functions_vol2.R")

#Sample size 
N <- 100

#Generate some survival data using one covariate beta

alfa <- 0.5
mu <- 0
sigma_g <- 2

#Effect sizes (arbitrary)
betas <- rnorm(200,0,sd=sqrt(sigma_g))
betas[sample(1:length(betas),size=190)] <- 0


maf <- runif(length(betas),0.4,0.6)
mafs <- c()
for(i in 1:length(betas)){
  mafs <- c(mafs,rep(maf[i],N))
}

#Generate genotypes
X <- matrix(rbinom(N*length(betas),2,mafs),ncol=length(betas))
X <- scale(X)

b <- (exp(mu+X%*%betas +rnorm(N,mean=0,sd=sqrt(0.2))))**(-1/alfa)  #Add also a small error

y <- rweibull(N,shape = alfa,scale=b)

#Add the intercept to the model matrix (corresponds to mu)
#X <- cbind(X)

#Prepare a data frame from those
data_set <- data.frame(y,X)

#Add randomly some failure times
data_set$research_end <- runif(N,25,30)

#Failure seen or not
data_set$failure <- (data_set$y < data_set$research_end)*1

data_set$time_in_research <- pmin(data_set$y,data_set$research_end)


#Construct a sampler doing sampling one by one from the distributions--------------------------------
#Fix the sample size
n <- 1000

#Prior values for alfa
alfa0 <- 6
kappa0 <- 1
#Prior values for sigma_g
sigma_g_alfa <- 0.001  #0.001 
sigma_g_beta <- 0.001  #0.001  #The prior is sigma_beta/2


#Prior values for beta
mu0 <- rep(0,length(betas))
sigma_g <- 1


#Set some starting values
BETA <- rep(0,length(betas))
alfa <- 5
#mu <- 0
prob <- 0.5   #Needs to be thought about

#BETA <- c(mu,BETA)  #First beta is the mu

#First gamma stands for intercept which is always 1, thus gamma_vec[1]=1 
gamma_vec <- sample(c(0,1),size=length(BETA),replace=T)
#gamma_vec <-  rep(1,length(BETA))

#Function for calculating the probability of including a SNP
gamma_i_prob_old <- function(prob,ind,g_vec,sigma_g){
  gamma_vec_0 <- copy(g_vec)
  gamma_vec_1 <- copy(g_vec)
  gamma_vec_0[ind] <- 0  #We shift it by one because the first (intercept) is always included
  gamma_vec_1[ind] <- 1
  
  D1 <- diag(data_set$failure)-diag(data_set$time_in_research**alfa)
  one_vec <- as.matrix(rep(1,N),ncol=1)
  
  X_0 <- (1/sigma_g) * diag(sum(gamma_vec_0)) + 2*t(X[,gamma_vec_0==1]) %*% diag(data_set$time_in_research**alfa) %*% X[,gamma_vec_0==1]
  X_1 <- (1/sigma_g) * diag(sum(gamma_vec_1)) + 2*t(X[,gamma_vec_1==1]) %*% diag(data_set$time_in_research**alfa) %*% X[,gamma_vec_1==1]
  
  S_gamma_0 <- X[,gamma_vec_0==1] %*% solve(X_0) %*% t(X[,gamma_vec_0==1])
  S_gamma_1 <- X[,gamma_vec_1==1] %*% solve(X_1) %*% t(X[,gamma_vec_1==1])
  
  #det_factor <- det(X_0)**(-0.5) * det(X_1)**(0.5)
  #Re-arrange X
  X_arranged <- X[,gamma_vec_0==1]
  X_arranged <- cbind(X_arranged,X[,ind])
  A <- (1/sigma_g)*diag(sum(gamma_vec_1)) + 2* t(X_arranged) %*% diag(data_set$time_in_research**alfa) %*% X_arranged
  det_factor <- sqrt(A[nrow(A),ncol(A)] - A[nrow(A),1:(ncol(A)-1)] %*% solve(A[1:(nrow(A)-1),1:(nrow(A)-1)]) %*% A[1:(ncol(A)-1),nrow(A)])
  
  return((1+ ((1-prob)/prob) * exp(t(one_vec) %*% D1 %*% (S_gamma_0-S_gamma_1) %*% D1 %*% one_vec) * det_factor) **(-1))
}

gamma_i_prob <- function(prob,ind,g_vec,sigma_g){
  #We shift the i to the last position of the matrices
  gamma_vec_0 <- copy(g_vec)
  gamma_vec_1 <- copy(g_vec)
  gamma_vec_0[ind] <- 0  #We shift it by one because the first (intercept) is always included
  gamma_vec_1[ind] <- 1
  
  D1 <- diag(data_set$failure)-diag(data_set$time_in_research**alfa)
  one_vec <- as.matrix(rep(1,N),ncol=1)
  
  #Re-arrange X so that the last column would be ith SNP
  X_arranged <- X[,gamma_vec_0==1]
  X_arranged <- cbind(X_arranged,X[,ind])
  
  X_0 <- (1/sigma_g) * diag(sum(gamma_vec_0)) + 2*t(X_arranged[,1:(ncol(X_arranged)-1)]) %*% diag(data_set$time_in_research**alfa) %*% X_arranged[,1:(ncol(X_arranged)-1)]
  X_1 <- (1/sigma_g) * diag(sum(gamma_vec_1)) + 2*t(X_arranged) %*% diag(data_set$time_in_research**alfa) %*% X_arranged
  
  X_0_inv <- solve(X_0)
  #C11 <- solve(X_0 - (X_1[ncol(X_1),ncol(X_1)])**(-1) * X_1[1:(nrow(X_1)-1),ncol(X_1)] %*% t(X_1[1:(nrow(X_1)-1),ncol(X_1)]))
  
  #Sherman-Morrison
  C11 <- X_0_inv + as.vector(X_1[nrow(X_1),ncol(X_1)]**(-1)/(1 - X_1[nrow(X_1),ncol(X_1)]**(-1) * t(X_1[1:(nrow(X_1)-1),ncol(X_1)] %*% X_0_inv %*% X_1[1:(nrow(X_1)-1),ncol(X_1)]))) *
    (X_0_inv %*% X_1[1:(nrow(X_1)-1),ncol(X_1)] %*% t(X_1[1:(nrow(X_1)-1),ncol(X_1)]) %*% X_0_inv)
  
  C22 <- (X_1[ncol(X_1),ncol(X_1)] - t(X_1[1:(nrow(X_1)-1),ncol(X_1)]) %*% X_0_inv %*% X_1[1:(nrow(X_1)-1),ncol(X_1)])**(-1)
  
  C21 <- -C22 %*% (X_1[1:(nrow(X_1)-1),ncol(X_1)]) %*% X_0_inv
  
  X_1_inv <- cbind(rbind(C11,C21),c(t(C21),C22))
  
  S_gamma_0 <- X_arranged[,1:(ncol(X_arranged)-1)] %*% X_0_inv %*% t(X_arranged[,1:(ncol(X_arranged)-1)])
  S_gamma_1 <- X_arranged %*% X_1_inv %*% t(X_arranged)
  
  #det_factor <- det(X_0)**(-0.5) * det(X_1)**(0.5)
  
  det_factor <- sqrt(X_1[nrow(X_1),ncol(X_1)] - X_1[nrow(X_1),1:(ncol(X_1)-1)] %*% X_0_inv %*% X_1[1:(ncol(X_1)-1),nrow(X_1)])
  
  return((1+ ((1-prob)/prob) * exp(t(one_vec) %*% D1 %*% (S_gamma_0-S_gamma_1) %*% D1 %*% one_vec) * det_factor) **(-1))
}


abscissae_set <- abscissae_setter_beta(alfa,BETA,sigma_g)

res_matrix <- matrix(rep(NA,(length(betas)+2)*n),ncol=length(betas)+2)
t1 <- Sys.time()
for(i in 1:n){
  #1) Update the shape parameter
  alfa <- ars(1,function(w) f_alfa_vec(w,BETA),function(w) f_alfa_derivative_vec(w,BETA),x=abscissae_set$alfa,m=5,lb=T,xlb=0)
  #plot(seq(3,7,0.1),f_alfa_vec(seq(3,7,0.1),LAMBDA),type="l")
  #alfa <- 5
  #3) Update gamma for each covariate (SNP marker)
  #If gamma = 1, then update beta  
  for(j in 1:length(betas)){
    gamma_prob <- gamma_i_prob(prob,j,gamma_vec,sigma_g)
  #  cat(paste0(j,": "))
  #  cat(gamma_prob)
  #  cat(" , ")
    gamma_sample <- rbinom(1,1,unlist(gamma_prob))
    gamma_vec[j] <- gamma_sample 
    if(gamma_sample==1){
      BETA[j] <- ars(1,function(w) f_beta_vec(w,alfa,ind=j,sigma_g),function(w) f_beta_derivative_vec(w,alfa,ind=j,sigma_g),x=abscissae_set$beta[,j],m=5)
    }
  }
  #4) Update the betas (including the intercept as the first beta) where gamma_vec = 1
  #Those where gamma is 0, BETA = 0
  BETA[gamma_vec==0] <- 0
  
  #Calculate new sigma_g value (using only the non 0)
  sigma_g <- rinvgamma(1,shape=sigma_g_alfa+0.5*sum(gamma_vec),rate=0.5 * (t(BETA[gamma_vec==1]) %*% BETA[gamma_vec==1]+sigma_g_beta))

  #Use abscissae finder only for some time
  if(i %% 10 ==0){
    t2 <- Sys.time()
    abscissae_set <- abscissae_setter_beta(alfa,BETA,sigma_g)
    #print(Sys.time()-t2)
  }
  if(i %% 5 ==0){
    print(i)
    print(Sys.time()-t1)
  }
  #print(paste0(i,": Alfa is ",new_alfa,". Beta is ",new_beta))
  res_matrix[i,] <- c(alfa,BETA,sigma_g)
  #res_matrix <- rbind(res_matrix,c(alfa,BETA,sigma_g))
}
print(Sys.time()-t1)
