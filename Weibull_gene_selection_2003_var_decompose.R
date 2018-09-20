rm(list=ls())
gc()

library(ars)
library(invgamma)
library(data.table)
library(MASS)

#Load the ARS functions
source("/Users/admin/Desktop/Surv1/ARS_functions_2003.R")



#Sample size 
N <- 1000

#Generate some survival data using one covariate beta

alfa <- 5
mu <- 0
sigma_g <- 1
sigma_e <- 0.01

#Effect sizes (arbitrary)
betas <- rnorm(30,0,sigma_g)
betas[sample(1:length(betas),size=25)] <- 0


maf <- runif(length(betas),0.4,0.6)
mafs <- c()
for(i in 1:length(betas)){
  mafs <- c(mafs,rep(maf[i],N))
}

#Generate genotypes
X <- matrix(rbinom(N*length(betas),2,mafs),ncol=length(betas))
X <- scale(X)

b <- (exp(mu+X%*%betas+rnorm(N,0,sd=sqrt(sigma_e))))**(-1/alfa)

y <- rweibull(N,shape = alfa,scale=b)

#Add the intercept to the model matrix (corresponds to mu)
X <- cbind(X)

#Prepare a data frame from those
data_set <- data.frame(y,X)

#Add randomly some failure times
data_set$research_end <- runif(N,25,30)

#Failure seen or not
data_set$failure <- (data_set$y < data_set$research_end)*1

data_set$time_in_research <- pmin(data_set$y,data_set$research_end)


#Construct a sampler doing sampling one by one from the distributions--------------------------------
#Fix the sample size
n <- 100

#Prior values for alfa
alfa0 <- 6
kappa0 <- 1
#Prior values for sigma_g
sigma_g_alfa <- 0.001  #0.001 
sigma_g_beta <- 0.001  #0.001  #The prior is sigma_beta/2

#Prior values for sigma_e
sigma_e_alfa <- 0.001  #0.001 
sigma_e_beta <- 0.001  #0.001  #The prior is sigma_beta/2


#Recommended value for variance factor
C <- 1 #100

#Prior values for beta
mu0 <- rep(0,length(betas))
sigma_g <- 1
sigma_e <- 0.1


#Set some starting values
BETA <- rep(0,length(betas))
alfa <- 5
#mu <- 0
prob <- 0.5   #Needs to be thought about

#BETA <- c(mu,BETA)  #First beta is the mu

#First gamma stands for intercept which is always 1, thus gamma_vec[1]=1 
gamma_vec <- rep(1,length(BETA))

#Function for calculating the probability of including a SNP
gamma_i_prob <- function(lambda,prob,ind,g_vec,sigma_e,sigma_g){
  gamma_vec_0 <- copy(g_vec)
  gamma_vec_1 <- copy(g_vec)
  gamma_vec_0[ind] <- 0  #We shift it by one because the first (intercept) is always included
  gamma_vec_1[ind] <- 1
  #S_gamma_0 <- sum(lambda**2) - (C/(1+C)) * t(lambda) %*% X[,gamma_vec_0==1 & c(F,rep(T,length(gamma_vec_0)-1))] %*% solve(t(X[,gamma_vec_0==1 & c(F,rep(T,length(gamma_vec_0)-1))]) %*% X[,gamma_vec_0==1 & c(F,rep(T,length(gamma_vec_0)-1))]) %*% t(X[,gamma_vec_0==1 & c(F,rep(T,length(gamma_vec_0)-1))]) %*% lambda
  #S_gamma_1 <- sum(lambda**2) - (C/(1+C)) * t(lambda) %*% X[,gamma_vec_1==1 & c(F,rep(T,length(gamma_vec_1)-1))] %*% solve(t(X[,gamma_vec_1==1 & c(F,rep(T,length(gamma_vec_1)-1))]) %*% X[,gamma_vec_1==1 & c(F,rep(T,length(gamma_vec_1)-1))]) %*% t(X[,gamma_vec_1==1 & c(F,rep(T,length(gamma_vec_1)-1))]) %*% lambda
  
  X_0 <- t(X[,gamma_vec_0==1]) %*% X[,gamma_vec_0==1]
  X_1 <- t(X[,gamma_vec_1==1]) %*% X[,gamma_vec_1==1]
  
  S_gamma_0 <- sum(lambda**2) - (C*sigma_g/(sigma_e+C*sigma_g)) * t(lambda) %*% X[,gamma_vec_0==1] %*% solve(X_0) %*% t(X[,gamma_vec_0==1]) %*% lambda
  S_gamma_1 <- sum(lambda**2) - (C*sigma_g/(sigma_e+C*sigma_g)) * t(lambda) %*% X[,gamma_vec_1==1] %*% solve(X_1) %*% t(X[,gamma_vec_1==1]) %*% lambda
  
  return((1+ ((1-prob)/prob) * sqrt(1/(C*sigma_g+sigma_e)) * exp(-(1/(2*sigma_e))*(S_gamma_0-S_gamma_1)) * det(X_0)**(-0.5) * det(X_1)**(0.5)) **(-1))
}

#Set the default value for lambda mu
#LAMBDA <- rep(mu,N)
LAMBDA <- X%*%betas

abscissae_set <- abscissae_setter(alfa,LAMBDA,0,sigma_e)

res_matrix <- matrix(rep(NA,(length(betas)+3)*n),ncol=length(betas)+3)
t1 <- Sys.time()
for(i in 1:n){
  #1) Update the shape parameter
  alfa <- ars(1,function(w) f_alfa_vec(w,LAMBDA),function(w) f_alfa_derivative_vec(w,LAMBDA),x=abscissae_set$alfa,m=5,lb=T,xlb=0)
  #plot(seq(3,7,0.1),f_alfa_vec(seq(3,7,0.1),LAMBDA),type="l")
  
  print(i)
  
  #2) Update the scale parameters per person
  for(j in 1:N){
    LAMBDA[j] <- ars(1,function(w) f_lambda_vec(w,j,alfa,sigma_e,gamma_vec,BETA),function(w) f_lambda_derivative_vec(w,j,alfa,sigma_e,gamma_vec,BETA),m=5,x=abscissae_set$lambda[j,])
  }
  
  #3) Update gamma for each covariate (SNP marker)
  for(j in 1:length(betas)){
    cat(gamma_i_prob(LAMBDA,prob,j,gamma_vec,sigma_e,sigma_g))
    cat(" , ")
    gamma_vec[j] <- rbinom(1,1,gamma_i_prob(LAMBDA,prob,j,gamma_vec,sigma_e,sigma_g))    #+1 because intercept is first
  }
  
  #4) Update the betas (including the intercept as the first beta) where gamma_vec = 1
  
  #Those were gamma is 0, BETA = 0
  BETA[gamma_vec==0] <- 0
  
  beta_mean <- sigma_g * C/(sigma_e+C*sigma_g) * solve(t(X[,gamma_vec==1])%*%X[,gamma_vec==1]) %*% t(X[,gamma_vec==1]) %*% LAMBDA
  beta_var <- (1/sigma_e+1/(C*sigma_g))**(-1)* solve(t(X[,gamma_vec==1])%*%X[,gamma_vec==1])
  
  if(F){
    sampled_values <- c()
    
    #We sample the last value and then condition on that
    while(length(beta_mean)>=2){
      #print(length(beta_mean))
      mu_2 <- beta_mean[length(beta_mean)]
      sigma_2 <- beta_var[length(beta_mean),length(beta_mean)]
      
      sampled_value <- rnorm(1,mean=mu_2,sd=sqrt(sigma_2))
      sampled_values <- c(sampled_value,sampled_values)
      #Renew mean and variance
      beta_mean <- beta_mean[1:(length(beta_mean)-1)] + beta_var[nrow(beta_var),1:(nrow(beta_var)-1)] * sigma_2 * (sampled_value-mu_2)
      
      beta_var <- beta_var[1:(nrow(beta_var)-1),1:(nrow(beta_var)-1)] - t(t(beta_var[1:(nrow(beta_var)-1),nrow(beta_var)])) %*% beta_var[nrow(beta_var),1:(nrow(beta_var)-1)] * (sigma_2)**(-1)
    }
    #Final sample
    sampled_values <- c(rnorm(1,mean=beta_mean,sd=sqrt(beta_var)),sampled_values)
    
    BETA[gamma_vec==1] <- sampled_values
  }
  BETA[gamma_vec==1] <- mvrnorm(1,beta_mean,beta_var)
  
  
  #Calculate new sigma_g value (using only the non 0)
  sigma_g <- rinvgamma(1,shape=sigma_g_alfa+0.5*sum(gamma_vec),rate=0.5 * ((1/C)*(t(BETA[gamma_vec==1]) %*% t(X[,gamma_vec==1]) %*% X[,gamma_vec==1] %*% BETA[gamma_vec==1])+sigma_g_beta))
  #Robust estimate:
#  sigma_g <- rinvgamma(1,shape=sigma_g_alfa+0.5*sum(gamma_vec),rate=0.5 * ((1/C)*(t(BETA[gamma_vec==1]) %*% BETA[gamma_vec==1])+sigma_g_beta))
#  sigma_g <- 1
  print(sigma_g)
  
  #Calculate new sigma_e value
  sigma_e <- rinvgamma(1,shape=sigma_e_alfa+0.5*N,rate=0.5 * (sum((LAMBDA-X[,gamma_vec==1]%*%BETA[gamma_vec==1])**2) +sigma_e_beta))
  print(sigma_e)
  #Use abscissae finder only for some time
  if(i %% 2 ==0){
    #print(i)
    t2 <- Sys.time()
    abscissae_set <- abscissae_setter(alfa,LAMBDA,1,sigma_e)
   # print(Sys.time()-t2)
  }
  
  #print(paste0(i,": Alfa is ",new_alfa,". Beta is ",new_beta))
  res_matrix[i,] <- c(alfa,BETA,sigma_g,sigma_e)
}
print(Sys.time()-t1)
