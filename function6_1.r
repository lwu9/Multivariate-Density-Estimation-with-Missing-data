setwd("\\\\wolftech.ad.ncsu.edu/cos/stat/Redirect/lwu9/Documents/740/FinalProject")
source("function5.R")
source("DataGenerate.R")
#############################
## The following does not use Gibbs, cannot run because of to large memory; can be used to verify Gibbs later
#library(gtools)
# posterior_sampler.lambda_sq <- function(data, knots, thetas) {
#   n=dim(data)[1]; p=dim(data)[2]; m=length(thetas)
#   # obtain the variance of each colunm in data, take it as one parameter of InvGamma prior 
#   a = rep(n**0.4+1, p); b = apply(data, 2, var, na.rm=T) 
#   # permutation with replacement to get qn's. The number of terms: m^n, 
#   # qns is a (m^n) by n matrix
#   qns = permutations(m,n, c(1:m), repeats.allowed=T)
#   # probabilities to determine which component to sample in the mixture Cauchy model
#   prob = sapply(1:m^n, function(i) prod(thetas[qns[i,]]))
#   K = sample.int(m^n, 1, prob=prob)
#   #knots: p by m matrix
#   for (i in 1:p) {
#     # new parameters for ith component in mixture Cauchy model
#     newa = n/2+a[i]+1
#     newb = sum((data[,i]-knots[i,qns[K,]])**2)+b[i]
#     # After fixing K, lambda_sq[i]'s are independent
#     lambda_sq[i] = 1/rgamma(1, newa, newb)
#   }
#   
#   return(lambda_sq)
# }

posterior_Gibbs_sampler.lambda_sq <- function(data, N.iter, N.burn, knots, thetas, a, b) {
  #knots: p by m matrix; thetas: length is m
  # n: sample size; p: dimension of data
  n = dim(data)[1]; p = dim(data)[2]
  m = length(thetas)
  #Compute total number of iterations to run:
  N.tot = N.iter + N.burn
  # set initial values
  lambda_sq = matrix(0, N.tot, p)
  lambda_sq[1,] = apply(data, 2, var, na.rm=TRUE)/n**(2/5)
  K = rep(1, n)
  
  # Start Gibbs sampling
  for (l in 2:N.tot) {
    lambda_sq.old = as.vector(lambda_sq[l-1,])
    for (j in 1:n) {
      prob = sapply(1:m, function(k) 
        thetas[k]*prod(dnorm(data[j,], knots[,k], sqrt(lambda_sq.old))) )
      # Sample K to determine which component in the jth sample
      K[j] = sample.int(m, 1, prob = prob)
    }
    newa = n/2 + a
    # length(K) is n, so knots[,K] is p by n matrix
    newb = apply((data - t(knots[,K]))**2, 2, sum)/2 + b
    lambda_sq.new=rep(1,p)
    for(i in 1:p){
        lambda_sq.new[i] = 1/rgamma(1, newa[i], newb[i])
    }

    lambda_sq[l,] = lambda_sq.new
  }
  
  #End of Gibbs sampling
  return(list(lambda.square=lambda_sq[(N.burn+1):N.tot,], knots=knots, theta=thetas))
}

#######################################
density_estimate_new <- function(data,m,N_burn=1000,init1,init2,N_iter,N_burn_out=0) {
  n=nrow(data)
  p=ncol(data)
  theta=matrix(0,N_iter+1,m)
  lambda_sq = matrix(0, N_iter+1, p)
  
  # initialize theta and lambda_sq
  theta[1,] = init1; lambda_sq[1,] = init2
  
  #record the position where the variable is NA
  #miss_ind: p by n matrix, TRUE means NA; 
  #miss_ind[i,j] = TRUE, i.e., data[j,i] = NA
  miss_ind = matrix(0, p, n)
  for (i in 1:p){
    miss_ind[i,] <- is.na(data[,i])
  }
  
  #Get a data set from "data" where both all variables (x, y,z,...) are observed
  com_data=data[complete.cases(data), ]
  com_n=nrow(com_data)
  
  #set parameters for prior of lambda_sq
  a = rep(com_n**0.4+1, p); b = apply(com_data, 2, var, na.rm=TRUE)
  
  # compute the knots and set the bandwidth for each variable
  # knots = matrix(0,p,m)
  # x=com_data[,1]
  # knots[1,] = c(min(x),sort(x)[round((1:(m-2))*com_n/(m-1))],max(x))
  # for (i in 2:p) {
  #  y <- com_data[,i]
  #  y.ord=y[order(x)]
  #  knots[i,]=y.ord[c(1,round((1:(m-2))*com_n/(m-1)),com_n)]
  # }

  knots = matrix(0,p,m);
  x=com_data[,1]
  index_order=round((1:(m-2))*com_n/(m-1))
  index_order=ifelse(index_order==0,1,index_order)
  knots[1,] = c(min(x),sort(x)[index_order],max(x));
  for (i in 2:p) {
    y <- com_data[,i]
    y.ord=y[order(x)]
    index_order=round((1:(m-2))*com_n/(m-1))
    index_order=ifelse(index_order==0,1,index_order)
    knots[i,]=y.ord[c(1,index_order,com_n)]

  }
  
  
  
  #OLD Method
  if(com_n==n){
    for(j in 1:N_iter){
      lambdas = sqrt(lambda_sq[j,])
      pGs.theta = posterior_Gibbs_sampler(data,m, 1, N_burn, knots, lambdas)
      theta[j+1,] = as.vector(pGs.theta[[1]])
      pGs.lambda_sq = posterior_Gibbs_sampler.lambda_sq(data, 1, N_burn, knots, theta[j+1,], a, b)
      lambda_sq[j+1,] = as.vector(pGs.lambda_sq[[1]])
      
    }  
    return(list(theta=theta[(N_burn_out+1):(N_iter+1),],knots=knots,
              lambda=sqrt(lambda_sq[(N_burn_out+1):(N_iter+1),])))
  }
  
  #NEW Method
  if(com_n!=n){
    for(j in 1:N_iter){
      lambdas = sqrt(lambda_sq[j,])
      for (i in 1:n){
        data[i,] = conditional_sampler(data[i,],miss_ind[,i],lambdas,knots,theta[j,])
      }
      pGs.theta = posterior_Gibbs_sampler(data,m, 1, N_burn, knots, lambdas)
      theta[j+1,] = as.vector(pGs.theta[[1]])
      pGs.lambda_sq = posterior_Gibbs_sampler.lambda_sq(data, 1, N_burn, knots, theta[j+1,], a, b)
      lambda_sq[j+1,] = as.vector(pGs.lambda_sq[[1]])
      # ps.lambda_sq = posterior_sampler.lambda_sq(data, knots, theta[j+1,])
      # lambda_sq[j+1,] = as.vector(ps.lambda_sq[[1]])
    }  
  }
  return(list(theta=theta[(N_burn_out+1):(N_iter+1),],knots=knots,
              lambda=sqrt(lambda_sq[(N_burn_out+1):(N_iter+1),])))
}

