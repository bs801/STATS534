#the inverse of the logit function
library(MASS)
inverseLogit <- function(x)
{
  return(exp(x)/(1+exp(x))); 
}

#function for the computation of the Hessian
inverseLogit2 <- function(x)
{
  return(exp(x)/(1+exp(x))^2); 
}

#computes pi_i = P(y_i = 1 | x_i)
getPi <- function(x,beta)
{
  x0 = cbind(rep(1,length(x)),x);
  
  return(inverseLogit(x0 %*% beta ));
}

#another function for the computation of the Hessian
getPi2 <- function(x,beta)
{
  x0 = cbind(rep(1,length(x)),x);
  return(inverseLogit2(x0 %*%beta));
}

#logistic log-likelihood (formula (3) in your handout)
logisticLoglik <- function(y,x,beta)
{
  Pi = getPi(x,beta);
  return(sum(y*log(Pi))+sum((1-y)*log(1-Pi)));
}

#hessian matrix for the posterior log likelihood l*
getHessian <- function(y,x,beta)
{
  hessian = matrix(0,2,2);
  Pi2 = getPi2(x,beta);
  
  hessian[1,1] = -1-sum(Pi2);
  hessian[1,2] = -sum(Pi2*x);
  hessian[2,1] = hessian[1,2];
  hessian[2,2] = -1-sum(Pi2*x^2);
  
  return(hessian);
}

#Gradient of posterior log likelihood l*
getGradient <- function(y,x,beta)
{
  gradient = matrix(0,2,1);
  Pi = getPi(x,beta);

  gradient[1,1] = -beta[1]+sum(y-Pi);
  gradient[2,1] = -beta[2]+sum((y-Pi)*x);
  
  return(gradient);
}

#Calculating the log determinant
logdet <- function(R){
  determinant = NULL;
  eigenvalues = eigen(R)$values;
  determinant = sum(log(eigenvalues))
  return(determinant);
}

#Newton Raphson method for finding the beta0
getcoefNR <- function(y,x,data)
{
  #2x1 matrix of coefficients`
  beta = matrix(0,2,1);
  epsilon = 1e-4
  
  #infinite loop unless we stop it someplace inside
  while( TRUE)
  {
    newBeta = beta - solve(getHessian(y,x,beta))%*%getGradient(y,x,beta);
    
    if(abs(beta[1]-newBeta[1]) < epsilon && abs(beta[2]-newBeta[2]) < epsilon){
      break;
    }
    beta = newBeta;

  }
  return(beta);
}

#This function computes the Laplace approximation (4) of the marginal likelihood
getLaplaceApprox  <- function(response, explanatory, data, betaMode){
    y = data[,response];
    x = data[,explanatory];
    
    log_PD = log(2*pi) + logisticLoglik(y,x,betaMode) - (1/2)* logdet((-1)*getHessian(y,x,betaMode))

    return(log_PD)
}

#The function returns Log likelihood of the posterior distribution
logisticLoglik2 <- function(y,x,betaMode)
{
  return(-log(2*pi)-(1/2)*sum(betaMode^2) + logisticLoglik(y,x,betaMode));
}

# The function calculate beta estimate from Metropolis-Hastings algorithm
getPosteriorMeans <- function(response,explanatory,data,betaMode,niter){
    y = data[,response];
    x = data[,explanatory];
    beta_list = list(c(betaMode[1,1],betaMode[2,1]))
    
    for(i in 1:niter){
        beta_prev = beta_list[[length(beta_list)]]
        beta_sample = MASS:::mvrnorm(n=1,mu=beta_prev,Sigma=(-1)*solve(getHessian(y,x,beta_prev)))

        
        LL_prev = logisticLoglik2(y,x,beta_prev)
        LL_sample = logisticLoglik2(y,x,beta_sample)
        if(LL_sample  > LL_prev  || LL_sample  == LL_prev){
            beta_list[[length(beta_list)+1]] = beta_sample
        }else{
            u <- runif(1)
            if(log(u) <  LL_sample-  LL_prev || log(u) ==  LL_sample-  LL_prev){
              beta_list[[length(beta_list)+1]] = beta_sample 
            }else{
              beta_list[[length(beta_list)+1]] = beta_prev
            }
        }
    }
    
    beta_list = beta_list[2:length(beta_list)];
    beta_matrix = matrix(unlist(beta_list),nrow=niter, byrow=TRUE)
    
    return(c(mean(beta_matrix[,1]), mean(beta_matrix[,2])))
}
# The function returns the MLE of beta0 and beta1
getcoefglm <- function(response,explanatory,data)
{
  return(coef(glm(data[,response] ~ data[,explanatory],family=binomial(link=logit))));
}


# main <- function(datafile){
#   
#   data = read.table(datafile,header = FALSE);
#   niter = 1e4
#   response = ncol(data);
#   explanatory = sample(1:response-1,1);
#   
#   betaMode = getcoefNR(data[,response],data[,explanatory],data);
#   betaMode = matrix(betaMode,nrow=2)
#   
#   q1 = getLaplaceApprox(response,explanatory,data,betaMode);
#   
#   q2 = getPosteriorMeans(response,explanatory,data,betaMode,niter);
#   
# 
#   beta_0_mean = mean(q2[,1])
#   beta_1_mean = mean(q2[,2])
#   
#   cat("The average beta_0 is ")
#   cat(beta_0_mean)
#   cat("\n")
#   cat("The average beta_1 is ")
#   cat(beta_1_mean)
#   cat("\n")
# }




#main("534binarydata.txt")


