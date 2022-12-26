bayesLogistic = function(apredictor,response,data,NumberOfIterations)
{
  #NOTE: YOU NEED THE PACKAGE 'RCDD' TO PROPERLY RUN THIS CODE
  #load the 'RCDD' package
  require(MASS);
  
  #these are your functions that need to be loaded before they can be called
  source('auxiliary_BeysianInference.R');
  
  #Calculate beta MLE
  betaMLE = getcoefglm(response,apredictor,data);
  betaMode = getcoefNR(data[,response],data[,apredictor],data);
  betaMode = matrix(betaMode,nrow=2);
  #Calculate log likelihood of posterior distribution
  logmarglik = getLaplaceApprox(response, apredictor, data, betaMode);
  #Calculate beta by Metropolis-Hastings algorithm
  betas = getPosteriorMeans(response,apredictor,data,betaMode,NumberOfIterations);
  
  
  return(list(apredictor=apredictor,logmarglik=logmarglik,beta0bayes=betas[1],beta1bayes=betas[2],beta0mle=betaMLE[1],beta1mle=betaMLE[2]))
}
#PARALLEL VERSION
#datafile = the name of the file with the data
#NumberOfIterations = number of iterations of the Metropolis-Hastings algorithm
#clusterSize = number of separate processes; each process performs one or more
#univariate regressions
main <- function(datafile,NumberOfIterations,clusterSize)
{
  #read the data
  data = read.table(datafile,header=FALSE);
  #the sample size is 148 (number of rows)
  #the explanatory variables are the first 60 columns for '534binarydata.txt'
  #the last column is the binary response
  response = ncol(data);
  lastPredictor = ncol(data)-1;
  #initialize a cluster for parallel computing
  cluster <- makeCluster(clusterSize, type = "SOCK")
  #run the bayesLogistic several times
  results = clusterApply(cluster, 1:lastPredictor, bayesLogistic,response,data,NumberOfIterations);
  # clusterExpand()
  #print out the results
  for(i in 1:lastPredictor)
  {
    cat('Regression of Y on explanatory variable ',results[[i]]$apredictor,
        ' has log marginal likelihood ',results[[i]]$logmarglik,
        ' with beta0 = ',results[[i]]$beta0bayes,' (',results[[i]]$beta0mle,')',
        ' and beta1 = ',results[[i]]$beta1bayes,' (',results[[i]]$beta1mle,')',
        '\n');
  }
  #destroy the cluster
  stopCluster(cluster);
}

require(parallel)

#main('534binarydata.txt',10000,10)