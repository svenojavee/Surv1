
#####################################################################

#Functions for Adaptive rejection sampling in Weibull survival context
#Date: 14.09.2018
#Author: Sven Erik Ojavee

#####################################################################


#Log alfa densities and derivatives-----------------------------------

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


#Log beta densities and derivatives-----------------------------------
#Function takes also the index (which beta are we sampling here)


f_beta <- function(beta,alfa,mu,ind,gamma_vec){
  BETA_USE <- copy(BETA)
  BETA_USE[ind] <- beta #Replace the specific beta index
  return(sum(data_set$failure*data_set[,1+ind]*beta-((data_set$time_in_research)**alfa)*exp(mu+as.matrix(X[,gamma_vec == 1]) %*% BETA_USE[gamma_vec==1]))-
           0.5*(1/sigma2)*(beta-mu0[ind])**2)
}
f_beta_vec <- function(beta,alfa,mu,ind,gamma_vec){
  return(unlist(lapply(beta,function(x) f_beta(x,alfa,mu,ind,gamma_vec))))
}

f_beta_derivative <- function(beta,alfa,mu,ind,gamma_vec){
  BETA_USE <- copy(BETA)
  BETA_USE[ind] <- beta #Replace the specific beta index
  return(sum(data_set$failure*data_set[,1+ind] - data_set$time_in_research**alfa*data_set[,1+ind] * exp(mu+as.matrix(X[,gamma_vec == 1]) %*% BETA_USE[gamma_vec == 1]))
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

f_beta_derivative_vec <- function(beta,alfa,mu,ind,gamma_vec){
  return(unlist(lapply(beta,function(x) f_beta_derivative(x,alfa,mu,ind,gamma_vec))))
}

#Log mu densities and derivatives-----------------------------------


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

#Log lambda densities and first derivative-----------------------------------
#Ind denotes a subject in the sample (not the SNP)

f_lambda <- function(lambda,ind,alfa,sigma2,gamma_vec,beta){
  return(sum(data_set$failure)*log(alfa) + data_set$failure[ind]*lambda+
           data_set$failure[ind]*(alfa-1)*log(data_set$time_in_research[ind])-
           exp(lambda)*data_set$time_in_research[i]**alfa +
           log(1/sqrt(sigma2)) -1/(2*sigma2)*(lambda-X[ind,gamma_vec] %*% beta[gamma_vec])**2)
}

f_lambda_derivative <- function(lambda,ind,alfa,sigma2,gamma_vec,beta){
  return(data_set$failure[ind]-
           exp(lambda)*data_set$time_in_research[i]**alfa +
           -1/(2*sigma2)*2*(lambda-X[ind,gamma_vec] %*% beta[gamma_vec]))
}

f_lambda_vec <- function(lambda,ind,alfa,sigma2,gamma_vec,beta){
  return(unlist(lapply(lambda,function(w) f_lambda(w,ind,alfa,sigma2,gamma_vec,beta) )))
}
f_lambda_derivative_vec <- function(lambda,ind,alfa,sigma2,gamma_vec,beta){
  return(unlist(lapply(lambda,function(w) f_lambda_derivative(w,ind,alfa,sigma2,gamma_vec,beta) )))
}



#Calculating the modes for ARS---------------------------------------

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
    beta_matrix[,i] <- beta_suggest[i] * c(1,1,1,1,1) + 10*cover2
  }
  alfa_suggest <- alfa_suggest * cover
  mu_suggest <- mu_suggest * c(1,1,1,1,1) + 10*cover2
  return(list(alfa=alfa_suggest,mu=mu_suggest,beta=beta_matrix))
}
