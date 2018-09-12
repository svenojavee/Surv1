rm(list=ls())
gc()

library(ars)
library(invgamma)
library(data.table)

#Sample size 
N <- 100

#Generate some survival data using one covariate beta

#Effect sizes (arbitrary)
betas <- rnorm(50,0,1)
betas2 <- rnorm(10,0,2)
alfa <- 5
mu <- -10

maf <- runif(length(betas),0.4,0.6)
mafs <- c()
for(i in 1:length(betas)){
  mafs <- c(mafs,rep(maf[i],N))
}

#Generate genotypes
X <- matrix(rbinom(N*length(betas),2,mafs),ncol=length(betas))
X <- scale(X)

#Generate scale parameters (in rweibull parametrisation)
#Multiplicative constant
K <- 1
b <- (exp(mu+X%*%betas))**(-1/alfa)*K

#Generate lifespans from Weibull with shape a=alfa and scale b
y <- rweibull(N,shape = alfa,scale=b)

#Prepare a data frame from those
data_set <- data.frame(y,X)

#Add randomly some failure times
data_set$research_end <- runif(N,5,20)

#Failure seen or not
data_set$failure <- (data_set$y < data_set$research_end)*1

data_set$time_in_research <- pmin(data_set$y,data_set$research_end)

#Log conditional posterior for alfa

#alfa has to be an element (not vector!)
f_alfa <- function(alfa,mu,beta){
  return((alfa0+sum(data_set$failure)-1)*log(alfa)+
           (alfa-1)*sum(data_set$failure*log(data_set$time_in_research))-
           sum((data_set$time_in_research)**alfa *exp(mu+as.matrix(data_set[,2:(length(betas)+1)]) %*% beta))
         -kappa0*alfa)
}

#alfa has to be an element (not vector!)
f_alfa_derivative <- function(alfa,mu,beta){
  return((alfa0+
            sum(data_set$failure)-1)/alfa +
           sum(data_set$failure*log(data_set$time_in_research)) -
           sum((data_set$time_in_research)**alfa*log(data_set$time_in_research) * exp(mu+as.matrix(data_set[,2:(length(betas)+1)]) %*% beta)) -
           kappa0)
}

#These are not vectorized versions
f_alfa_derivative_local <- function(alfa,mu,beta_use){
  return((alfa0+
            sum(data_set$failure)-1)/alfa +
           sum(data_set$failure*log(data_set$time_in_research)) -
           sum((data_set$time_in_research)**alfa*log(data_set$time_in_research) * exp(mu+as.matrix(data_set[,2:(length(betas)+1)]) %*% beta_use)) -
           kappa0)
}

f_alfa_derivative2 <- function(alfa,mu,beta_use){
  return(-(alfa0+sum(data_set$failure)-1)/(alfa**2)-
           sum((data_set$time_in_research)**alfa * log(data_set$time_in_research)**2 * exp(mu+as.matrix(data_set[,2:(length(betas)+1)]) %*% beta_use)))
}


#These are vectorized versions
f_alfa_vec <- function(alfa,mu,beta){
  return(unlist(lapply(alfa,function(x) f_alfa(x,mu,beta) )))
}
f_alfa_derivative_vec <- function(alfa,mu,beta){
  return(unlist(lapply(alfa,function(x) f_alfa_derivative(x,mu,beta) )))
}


###############
#Function takes also the index (which beta are we sampling here)
f_beta <- function(beta,alfa,mu,ind){
  BETA_USE <- copy(BETA)
  BETA_USE[ind] <- beta #Replace the specific beta index
  return(sum(data_set$failure*data_set[,1+ind]*beta-((data_set$time_in_research)**alfa)*exp(mu+as.matrix(data_set[,2:(length(betas)+1)]) %*% BETA_USE))-
           0.5*(1/sigma2)*(beta-mu0[ind])**2)
}
f_beta_vec <- function(beta,alfa,mu,ind){
  return(unlist(lapply(beta,function(x) f_beta(x,alfa,mu,ind))))
}

f_beta_derivative <- function(beta,alfa,mu,ind){
  BETA_USE <- copy(BETA)
  BETA_USE[ind] <- beta #Replace the specific beta index
  return(sum(data_set$failure*data_set[,1+ind] - data_set$time_in_research**alfa*data_set[,1+ind] * exp(mu+as.matrix(data_set[,2:(length(betas)+1)]) %*% BETA_USE))
         -1/(sigma2)*(beta-mu0[ind]))
}

f_beta_derivative_local <- function(beta,beta_use,alfa_use,mu_use,ind){
  beta_use[ind] <- beta #Replace the specific beta index
  return(sum(data_set$failure*data_set[,1+ind] - data_set$time_in_research**alfa_use*data_set[,1+ind] * exp(mu_use+as.matrix(data_set[,2:(length(betas)+1)]) %*% beta_use))
         -1/(sigma2)*(beta-mu0[ind]))
}
f_beta_derivative2 <- function(beta,beta_use,alfa_use,mu_use,ind){
  beta_use[ind] <- beta #Replace the specific beta index
  return(-sum(data_set$y**alfa_use * exp(mu_use+as.matrix(data_set[,2:(length(betas)+1)]) %*% beta_use) * (alfa_use*data_set[,1+ind])**2)-(1/sigma2))
}

f_beta_derivative_vec <- function(beta,alfa,mu,ind){
  return(unlist(lapply(beta,function(x) f_beta_derivative(x,alfa,mu,ind))))
}

################

f_mu <- function(mu,alfa_use,beta_use){
  return(mu*sum(data_set$failure)-sum(exp(mu+as.matrix(data_set[,2:(length(betas)+1)]) %*% beta_use)*data_set$time_in_research**alfa_use)-(1/(2*sigma_mu))*mu**2)
}
f_mu_derivative <- function(mu,alfa_use,beta_use){
  return(sum(data_set$failure)-sum(exp(mu+as.matrix(data_set[,2:(length(betas)+1)]) %*% beta_use)*data_set$time_in_research**alfa_use)-(1/(sigma_mu))*mu)
}
f_mu_derivative2 <- function(mu,alfa_use,beta_use){
  return(sum(exp(mu+as.matrix(data_set[,2:(length(betas)+1)]) %*% beta_use)*data_set$time_in_research**alfa_use)-(1/(sigma_mu)))
}
f_mu_vec <- function(mu,alfa_use,beta_use){
  return(unlist(lapply(mu,function(x) f_mu(x,alfa_use,beta_use))))
}
f_mu_derivative_vec <- function(mu,alfa_use,beta_use){
  return(unlist(lapply(mu,function(x) f_mu_derivative(x,alfa_use,beta_use))))
}

alfa_mode <- function(init_val,mu_use,beta_use,error){
  xi <- init_val
  xi_1 <- init_val+error+0.001
  while(abs(xi-xi_1)>error){
    xi_1 <- xi
    xi <- xi_1 - f_alfa_derivative_local(xi_1,mu_use,beta_use)/f_alfa_derivative2(xi_1,mu_use,beta_use)
    if(xi<0){
      return(init_val) #Give up
    }
  }
  return(xi)
}

mu_mode <- function(init_val,alfa_use,beta_use,error){
  init_val <- mu
  xi <- init_val
  xi_1 <- init_val+error+0.001
  counter <- 0
  
  while(abs(xi-xi_1)>error){
    counter <- counter + 1
    if(counter > 11) return(init_val)  #Failure
    xi_1 <- xi
    xi <- xi_1 - f_mu_derivative(xi_1,alfa_use,beta_use)/f_mu_derivative2(xi_1,alfa_use,beta_use)
    if(is.na(xi)) return(init_val)
  }
  return(xi)
}

beta_mode <- function(init_val,ind,alfa_use,mu_use,beta_use,error){
  xi <- init_val
  xi_1 <- init_val+error+0.001
  counter <- 0
  
  while(abs(xi-xi_1)>error){
    counter <- counter + 1
    if(counter > 11) return(init_val)  #Failure
    xi_1 <- xi
    xi <- xi_1 - f_beta_derivative_local(xi_1,beta_use,alfa_use,mu_use,ind)/f_beta_derivative2(xi_1,beta_use,alfa_use,mu_use,ind)
  }
  return(xi)
}

abscissae_setter <- function(alfa,beta,mu,sigma2){
  cover <- c(0.25,0.5,1,2,4)
  cover2 <- c(-2,-1,0,1,2)
  
  alfa_suggest <- alfa_mode(alfa,mu,beta,0.001) 
  
  mu_suggest <- mu_mode(mu,alfa,beta,0.001)
  
  beta_suggest <- copy(beta)
  beta_matrix <- matrix(rep(NA,length(beta)*length(cover)),nrow=5)
  for(i in 1:length(beta_suggest)){
    beta_suggest[i] <- beta_mode(beta_suggest[i],i,alfa_suggest,mu_suggest,beta_suggest,0.001)
    beta_matrix[,i] <- beta_suggest[i] * c(1,1,1,1,1) + 10*cover2*sigma2
  }
  alfa_suggest <- alfa_suggest * cover
  mu_suggest <- mu_suggest * c(1,1,1,1,1) + 10*cover2*sigma2
  return(list(alfa=alfa_suggest,mu=mu_suggest,beta=beta_matrix))
}


#The functions are really log-concave!
#plot(seq(0,25,0.1),f_alfa_vec(seq(0,25,0.1),beta=c(0,0,0)),type="l")
#BETA <-c(2,3,2)
#plot(seq(-10,-2,0.1),f_beta_vec(seq(-10,-2,0.1),1),type="l")
#plot(seq(-2,2,0.1),f_beta_vec(seq(-2,2,0.1),2),type="l")


#Construct a sampler doing sampling one by one from the distributions--------------------------------
#Fix the sample size
n <- 200

#Prior values for alfa
alfa0 <- 6
kappa0 <- 1
#Prior values for sigma2
sigma_alfa <- 0.001  #0.001
sigma_beta <- 0.001  #0.001

#Prior values for beta
mu0 <- rep(0,length(betas))
sigma2 <- 1

#Prior variance for mu
sigma_mu <- 10 #Just fix some big value


#Set some starting values
BETA <- rep(0,length(betas))
alfa <- 3
mu <- 0

res_matrix <- matrix(rep(NA,(length(betas)+3)*n),ncol=length(betas)+3)
abscissae_set <- abscissae_setter(alfa,BETA,mu,sigma2)
t1 <- Sys.time()
for(i in 1:n){
  alfa <- ars(1,function(w) f_alfa_vec(w,mu,BETA),function(w) f_alfa_derivative_vec(w,mu,BETA),x=abscissae_set$alfa,m=5,lb=T,xlb=0)
  print(i)
  mu <- ars(1,function(w) f_mu_vec(w,alfa,BETA),function(w) f_mu_derivative_vec(w,alfa,BETA),x=abscissae_set$mu,m=5)
  for(ind in 1:length(BETA)){
    #abscissae_set <- abscissae_setter(alfa,BETA,z)
    #print(ind)
    new_beta <- ars(1,function(w) f_beta_vec(w,alfa,mu,ind),function(w) f_beta_derivative_vec(w,alfa,mu,ind),x=abscissae_set$beta[,ind],m=5)
    BETA[ind] <- new_beta
    
  }
  #Calculate new sigma2 value
  sigma2 <- rinvgamma(1,shape=sigma_alfa+0.5*length(betas),scale=1/(0.5* sum(BETA**2)+sigma_beta))  
  
  #Use abscissae finder only for some time
  if(i<20){
    #print(i)
    t2 <- Sys.time()
    abscissae_set <- abscissae_setter(alfa,BETA,mu,sigma2)
    print(Sys.time()-t2)
  }
  
  #print(paste0(i,": Alfa is ",new_alfa,". Beta is ",new_beta))
  res_matrix[i,] <- c(alfa,mu,BETA,sigma2)
}


