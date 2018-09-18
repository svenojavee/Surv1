
#####################################################################

#Functions for Adaptive rejection sampling in Weibull survival context
#Date: 14.09.2018
#Author: Sven Erik Ojavee

#####################################################################


#Log alfa densities and derivatives-----------------------------------

#alfa has to be an element (not vector!)
f_alfa <- function(alfa,lambda_vec){
  return((alfa0+sum(data_set$failure)-1)*log(alfa)+
           (alfa-1)*sum(data_set$failure*log(data_set$time_in_research))-
           sum((data_set$time_in_research)**alfa *exp(lambda_vec))
         -kappa0*alfa)
}


#alfa has to be an element (not vector!)
f_alfa_derivative <- function(alfa,lambda_vec){
  return((alfa0+
            sum(data_set$failure)-1)/alfa +
           sum(data_set$failure*log(data_set$time_in_research)) -
           sum((data_set$time_in_research)**alfa*log(data_set$time_in_research) * exp(lambda_vec)) -
           kappa0)
}

#These are not vectorized versions
f_alfa_derivative_local <- function(alfa,lambda_vec){
  return((alfa0+
            sum(data_set$failure)-1)/alfa +
           sum(data_set$failure*log(data_set$time_in_research)) -
           sum((data_set$time_in_research)**alfa*log(data_set$time_in_research) * exp(lambda_vec)) -
           kappa0)
}

f_alfa_derivative2 <- function(alfa,lambda_vec){
  return(-(alfa0+sum(data_set$failure)-1)/(alfa**2)-
           sum((data_set$time_in_research)**alfa * log(data_set$time_in_research)**2 * exp(lambda_vec)))
}


#These are vectorized versions
f_alfa_vec <- function(alfa,lambda_vec){
  return(unlist(lapply(alfa,function(x) f_alfa(x,lambda_vec) )))
}
f_alfa_derivative_vec <- function(alfa,lambda_vec){
  return(unlist(lapply(alfa,function(x) f_alfa_derivative(x,lambda_vec) )))
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





#Log lambda densities and first derivative-----------------------------------
#Ind denotes a subject in the sample (not the SNP)

f_lambda <- function(lambda,ind,alfa,sigma2,gamma_vec,beta){
  return(data_set$failure[ind]*lambda-
           exp(lambda)*data_set$time_in_research[ind]**alfa  -1/(2*sigma2)*(lambda-X[ind,gamma_vec==1] %*% BETA[gamma_vec==1])**2)
}

f_lambda_derivative <- function(lambda,ind,alfa,sigma2,gamma_vec,beta){
  return(data_set$failure[ind]-
           exp(lambda)*data_set$time_in_research[ind]**alfa 
           -1/(2*sigma2)*2*(lambda-X[ind,gamma_vec==1] %*% beta[gamma_vec==1]))
}

f_lambda_vec <- function(lambda,ind,alfa,sigma2,gamma_vec,beta){
  return(unlist(lapply(lambda,function(w) f_lambda(w,ind,alfa,sigma2,gamma_vec,beta) )))
}
f_lambda_derivative_vec <- function(lambda,ind,alfa,sigma2,gamma_vec,beta){
  return(unlist(lapply(lambda,function(w) f_lambda_derivative(w,ind,alfa,sigma2,gamma_vec,beta) )))
}



#Calculating the modes for ARS---------------------------------------

alfa_mode <- function(init_val,lambda_use,error){
  xi <- init_val
  xi_1 <- init_val+error+0.001
  while(abs(xi-xi_1)>error){
    xi_1 <- xi
    xi <- xi_1 - f_alfa_derivative_local(xi_1,lambda_use)/f_alfa_derivative2(xi_1,lambda_use)
    if(xi<0){
      return(init_val) #Give up
    }
  }
  return(xi)
}


abscissae_setter <- function(alfa,lambda,mu,sigma2){
  cover <- c(0.25,0.5,1,2,4)
  cover2 <- c(-2,-1,0,1,2)
  
  alfa_suggest <- alfa_mode(alfa,lambda,0.001) 
  
  #Lambdas are heavily dependant on the mu values. Let's use this knowledge
  mu_temp <- rep(mu,length(lambda))
  lambda_matrix <- cbind(mu_temp-20*sigma2,mu_temp-5*sigma2,mu_temp,mu_temp+5*sigma2,mu_temp+10*sigma2)
  
  alfa_suggest <- alfa_suggest * cover
  return(list(alfa=alfa_suggest,lambda=lambda_matrix))
}
