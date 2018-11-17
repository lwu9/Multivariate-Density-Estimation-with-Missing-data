##########Functions
#Generate a r.v. from Dirichlet distribution
rdirichlet=function(n=1,alpha){
  if(min(alpha)<0){return("alpha must be positive")}
  m <- length(alpha)
  x <- matrix(rgamma(m*n,alpha),ncol=m,byrow=T)
  return(x/apply(x,1,sum))
}

####################
posterior_Gibbs_sampler<-function(data, m, N.iter,N.burn, knots, lambdas){
  n=dim(data)[1]
  #p: the dimension of variables
  p = length(lambdas)
  #set the parameter of Dirichlet prior:
  a=rep(1,m)
  #compute the knots # the middle numbers are equally-spaced ??
  #x.knot=c(min(x),sort(x)[round((1:(m-2))*n/(m-1))],max(x))
  #y.ord=y[order(x)]
  #y.knot=y.ord[c(1,round((1:(m-2))*n/(m-1)),n)]
  
  #set the bandwiths for x and y data:
  #lambda.x=sd(x)/n^(1/5); lambda.y=sd(y)/n^(1/5) ??
  #lambda.x=density(x)$bw; lambda.y=density(y)$bw
  #compute the joint density values at data and knots
  #knots: p by m matrix; lambdas: length is p
  f=matrix(0,n,m)
  for(k in 1:m){
    f[,k] = rep(1,n)
    for (j in 1:p){
      f[,k]=f[,k]*dnorm(data[,j],mean=knots[j,k],sd=lambdas[j])
    }
  }
  
  #Compute total number of iterations to run:
  N.tot=N.iter+N.burn
  #set initial values:
  theta=matrix(0,N.tot,m); K=1; theta[1,]=rep(1,m)/m
  
  #Start Gibbs sampling: # posterier distribution ??
  for(l in 2:N.tot){
    theta.old=as.vector(theta[l-1,])
    for(i in 1:n){
      fvalue=as.vector(f[i,])
      K[i]=sample.int(m,1,prob=theta.old*fvalue)
    }
    theta.new=rdirichlet(1,a+tabulate(K,nbins=m))
    theta[l,]=theta.new
  }
  #End of Gibbs sampling
  
  return(list(theta[(N.burn+1):N.tot,],knots, lambdas))
}

################### #
#sample y_miss (x_miss) given theta, x_obs, y_obs and x_miss (y_miss)

conditional_sampler <- function(vec, missx_ind, lambdas, knots, theta1) {
  p=length(lambdas)
  ind = c(1:p)[missx_ind==TRUE]; ind_Nomiss = c(1:p)[missx_ind==FALSE]
  if(length(vec)==length(ind_Nomiss)){return(vec)}
  a=NULL;
  m=ncol(knots)
  for (i in 1:m) {
    nomis=1
    for (l in ind_Nomiss){nomis=nomis*dnorm(vec[l],knots[l,i],lambdas[l])}
    a=c(a, theta1[i]*nomis)
  }
  beta=a/sum(a)
  k=sample(1:m,1,prob=beta)
  for (j in ind) {
    vec[j] <- rnorm(1,knots[j,k],lambdas[j])
  }
  return(vec)
}
# conditional_sampler<-function(x,lambda.x,lambda.y,x.knot,y.knot,theta){
#   a=NULL
#   m=length(x.knot)
#   for(i in 1:length(theta)){
#     a=c(a,theta[i]*dnorm(x,x.knot[i],lambda.x))
#   }
#   beta=a/sum(a)
#   k=sample(1:m,1,prob=beta)
#   return(rnorm(1,y.knot[k],lambda.y))
# }

############
density_estimate<-function(data,m,N_burn=1000,init,N_iter,N_burn_out=0){
  theta=matrix(0,N_iter+1,m)
  #x.knot=matrix(0,N_iter+1,m)
  #y.knot=matrix(0,N_iter+1,m)
  #lambda.x=matrix(0,N_iter+1,1)
  #lambda.y=matrix(0,N_iter+1,1)
  
  #initialize theta
  theta[1,]=init
  n=nrow(data)
  p=ncol(data)
  
  #record the position where x or y is NA
  #miss_ind: p by n matrix, TRUE means NA
  miss_ind = matrix(0, p, n)
  for (i in 1:p){
    miss_ind[i,] <- is.na(data[,i])
  }
  
  #Get a data set from "data" where both x and y are observed
  com_data=data[complete.cases(data), ]
  com_n=nrow(com_data)
  
  n=com_n 
  
  #compute the knots and set the bandwidth for each variable
  knots = matrix(0,p,m); lambdas = rep(0, p)
  x=com_data[,1]
  knots[1,] = c(min(x),sort(x)[round((1:(m-2))*n/(m-1))],max(x)); lambdas[1] = density(x)$bw
  for (i in 2:p) {
    y <- com_data[,i]
    y.ord=y[order(x)]
    knots[i,]=y.ord[c(1,round((1:(m-2))*n/(m-1)),n)]
    lambdas[i] = density(y)$bw
  }
  
  #OLD METHOD
  if(nrow(com_data)==nrow(data)){
    pGs=posterior_Gibbs_sampler(data, m, N_iter-N_burn_out, N_burn_out, knots, lambdas)
    return(list(theta=pGs[[1]],knots=knots,lambdas=lambdas))
  }
  #initialize x.knot and y.knot
  #x.knot[1,]=c(min(com_data[,1]),sort(com_data[,1])[round((1:(m-2))*com_n/(m-1))],max(com_data[,1]))
  #y.ord=com_data[,2][order(com_data[,1])]
  #y.knot[1,]=y.ord[c(1,round((1:(m-2))*com_n/(m-1)),com_n)]
  
  #initialize lambda.x and lambda.y
  #lambda.x[1,]=density(com_data[,1])$bw
  #lambda.y[1,]=density(com_data[,2])$bw
  
  #NEW METHOD
  if(nrow(com_data)!=nrow(data)){
    for(j in 1:N_iter){
      for(i in 1:nrow(data)){
        #sample missing values from the conditional distribution
        data[i,]=conditional_sampler(data[i,],miss_ind[,i],lambdas,knots,theta[j,])
        # if(a[i]==1){
        #   data[i,1]=conditional_sampler(data[i,2],lambda.y,lambda.x,y.knot,x.knot,theta[j,])
        #   missing_x[j,l1]=data[i,1]
        #   l1=l1+1
        # }
        # if(b[i]==1){
        #   data[i,2]=conditional_sampler(data[i,1],lambda.x,lambda.y,x.knot,y.knot,theta[j,])
        #   missing_y[j,l2]=data[i,2]
        #   l2=l2+1
        # }
      }
      
      #sample theta from the posterior distribution
      pGs=posterior_Gibbs_sampler(data,m,1,N_burn, knots, lambdas)
      theta[j+1,]=as.vector(pGs[[1]])
      #x.knot[j+1,]=c[[2]]
      #y.knot[j+1,]=c[[3]]
      #lambda.x[j+1,]=c[[4]]
      #lambda.y[j+1,]=c[[5]]
      
    }
  }
  return(list(theta=theta[(N_burn_out+1):(N_iter+1),],knots=knots,lambdas=lambdas))
}

## Some Summary functions
#marginal density of x
single_den=function(x,theta,x.knot,lambda.x){
  return(sum(theta*dnorm(x,mean=x.knot,sd=lambda.x)))
}

# vec: the row data, x is jth col, to get the condi expectation of x given the other variables
mu <- function(vec, j, theta, knots, lambdas) {
  ind = c(1:length(vec))[is.na(vec)]
  nomiss = c(1:length(vec))[-ind]
  den_prod = 1
  for (q in nomiss){
    den_prod = den_prod * dnorm(vec[q], knots[q,], lambdas[q])
  }
  weight = theta*den_prod
  return(sum(knots[j,]*weight)/sum(weight))
}


summaries <- function(data, j, theta, knots, lambdas) {
  x=data[,j]
  x.grid=seq(min(x),max(x),l=length(x));
  fx.hat=matrix(0,length(x.grid),nrow(theta))
  #mu.hat=matrix(0,length(x.grid),nrow(theta))
  for(i in 1:length(x.grid)){
    for(l in 1:nrow(theta)){
      #mu.hat[i,l]=mu(c(x.grid[i],NA),2,theta[l,],knots,lambdas)
      fx.hat[i,l]=single_den(x.grid[i],theta[l,],knots[j,],lambdas[j])
      #sigma.hat[i,l]=sqrt(lambda.y^2+mu2(x.grid[i],theta[l,],x.knot,y.knot,lambda.x)-mu.hat[i,l]^2)
    }
  }
  fx.hat.gibbs.l=apply(fx.hat,1,quantile,prob=0.025)
  fx.hat.gibbs=apply(fx.hat,1,mean)
  fx.hat.gibbs.u=apply(fx.hat,1,quantile,prob=0.975)
  return(list(fx.hat=fx.hat,fx.hat.gibbs.l=fx.hat.gibbs.l,fx.hat.gibbs=fx.hat.gibbs,fx.hat.gibbs.u=fx.hat.gibbs.u,x.grid=x.grid))
}

#Calculate the residuals of the prediction and true data using the new method
#Prediction of y
pred <- function(aa, j, theta, knots, lambdas) {
  data_miss = aa[[2]]
  # the missing index of the variable Xj
  xj_miss_ind = c(1:nrow(data_miss))[is.na(data_miss[,j])]
  if(length(xj_miss_ind)==0){
    return(list(true=NULL,xj.pred=NULL))
  }
  xj.pred=matrix(0,length(xj_miss_ind),nrow(theta))
  for(i in 1:length(xj_miss_ind)){
    for(l in 1:nrow(theta)){
      xj.pred[i,l]=mu(data_miss[xj_miss_ind[i],], j, theta[l,], knots, lambdas)
    }
  }
  return(list(true=aa[[1]][,j][xj_miss_ind], xj.pred=xj.pred))
}


errors <- function(aa, theta, knots, lambdas) {
  error1 = NULL
  for (j in 1:length(lambdas)) {
    pred_true = pred(aa,j,theta,knots, lambdas)
    if(length(pred_true[[1]])==0){
      error1=c(error1,0)
    }
    else{
      y.pred = pred_true$xj.pred
      y.pred.gibbs.l=apply(y.pred,1,quantile,prob=0.025)
      y.pred.gibbs=apply(y.pred,1,mean)
      y.pred.gibbs.u=apply(y.pred,1,quantile,prob=0.975)
      true = pred_true$true
      error1 = c(error1, sum((true-y.pred.gibbs)^2))
    }
  }
  return(error1)
}



cv_miss=function(data,index_pred){
  #the length of true data
  n=dim(data)[1]
  #the dim of parameter
  p=dim(data)[2]
  #the index of missing data
  index1=which(apply(data,1,function(x) any(is.na(x))))
  #the length of missing data
  n1=length(index1)
  #the length of complete data
  n2=n-n1
  # the test data set contains all comlete data
  
  if(n2<10)
    stop("the size of the complete data is smaller than 10")
  #size of the test test
  
  size=floor(n2/10)
  #size of the complete traning data
  size_train=n2-size
  
  #the index of complete data
  index=1:n
  if(n2==n)
    index2=1:n
  else
    index2=index[!index%in%index1]
  #split the complete data into k group
  group_index=list(10)
  for(i in 1:10){
    group_index[[i]]=index2[((i-1)*size+1):(i*size)]
  }
  
  squared_error=function(index,m){
    test_data1=data[index,]
    if(size==1)
      test_data1=matrix(test_data1,ncol=p)
    train_data=data[-index,]
    b=density_estimate(train_data,m,N_burn=100,init=rep(1/m,m),N_iter=100)
    knots=b$knots
    theta=b$theta
    lambdas=b$lambdas
    #leave one out predict
    predict_1=function(k){
      if(p==2){
        knots1=matrix(knots[-k,],1)
        lambdas1=lambdas[-k]
        test_data2=matrix(test_data1[,-k],length(test_data1[,-k]))
      }
      else{
        knots1=knots[-k,]
        
        lambdas1=lambdas[-k]
        test_data2=test_data1[,-k]
        if(nrow(test_data1)==1)
          test_data2=matrix(test_data2,ncol=p-1)
      }
      l=nrow(test_data2)
      predict_value=rep(0,nrow(test_data1))
      for(i in 1:nrow(test_data1)){
        
        weight1=apply(knots1,2,function(x) prod(dnorm(test_data2[i,],mean=x,sd=lambdas1)))
        value=apply(theta,1,function(x) sum(weight1*x/sum(weight1*x)*knots[k,]))
        predict_value[i]=mean(value)
      }
      return(mean((test_data1[,k]-predict_value)^2))
    }
    error=unlist(lapply(index_pred,predict_1))
    return(mean(error))
  }
  mean_square_error=function(m){
    error1=unlist(lapply(group_index,function(x) squared_error(x,m)))
    return(mean(error1))
  }
  #choose the possible value of m
  center=ceiling(size_train/log(size_train))
  seq1=(center-2):(center+2)
  step=floor((size_train-center-2)/3)
  seq3=seq(center+2+step,center+2+step*3,step)
  step=floor((center-3-2)/2)
  seq2=seq(center-2-2*step,center-2-step,step)
  m1=c(seq2,seq1,seq3)
  
  
  error2=unlist(lapply(m1,mean_square_error))
  return(m1[which.min(error2)])
  
}