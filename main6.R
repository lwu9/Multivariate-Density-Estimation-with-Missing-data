#Generate a raw data set and a missing data set
set.seed(740)
#set true mean and sd functions
mu.true=function(x){exp(x/6)-x+log(x^4+1)}
sigma.true=function(x){(x^2)*exp(-abs(x))}
#par(mfrow=c(2,2))
#curve(mu.true,-5,5); curve(sigma.true,-5,5)

#generate data:
n=100; x=rnorm(n,0,2)
y=mu.true(x)+sigma.true(x)*rnorm(n)
miss_x_prop = 0.4; miss_y_prop= 0.4
miss_x_ind = sample(1:n, size =round(miss_x_prop*n),  replace = FALSE)
miss_y_ind = sample((1:n)[-miss_x_ind], size =round(miss_y_prop*n),  replace = FALSE)
x_miss = x; y_miss = y
x_miss[miss_x_ind] = NA; y_miss[miss_y_ind] = NA
aa = list(cbind(x,y), cbind(x_miss, y_miss))

#Estimate theta from the missing data set
n=nrow(aa[[2]][complete.cases(aa[[2]]),])
m=ceiling(n/log(n))
m=10
b=density_estimate(aa[[2]],m,100,rep(1/m,m),100,0)

#########NEW METHOD
theta=b[[1]]; knots=b[[2]]; lambdas=b[[3]];

#########NEW METHOD WITH LAMBDA
com_data=aa[[2]][complete.cases(aa[[2]]), ]
p=ncol(aa[[2]])
init2=rep(1,p)
init1 = rep(1/m,m);
for (i in 1:p) {
  y <- com_data[,i]
  init2[i] = density(y)$bw^2
}

c=density_estimate_new(aa[[2]],m,N_burn=100,init1,init2,N_iter=100,N_burn_out=0)
theta_new=c[[1]]; knots_new=c[[2]]; lambdas_new=c[[3]]


#OLD METHOD
#####################
N.iter=100; N.burn=0
data=aa[[2]][complete.cases(aa[[2]]), ]
old_theta=posterior_Gibbs_sampler(data, m, N.iter,N.burn, knots, lambdas)[[1]]

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

summaries_new <- function(data, j, theta, knots, lambdas_new) {
  x=data[,j]
  x.grid=seq(min(x),max(x),l=length(x));
  fx.hat=matrix(0,length(x.grid),nrow(theta))
  #mu.hat=matrix(0,length(x.grid),nrow(theta))
  for(i in 1:length(x.grid)){
    for(l in 1:nrow(theta)){
      #mu.hat[i,l]=mu(c(x.grid[i],NA),2,theta[l,],knots,lambdas)
      fx.hat[i,l]=single_den(x.grid[i],theta[l,],knots[j,],lambdas_new[l,j])
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

pred_new <- function(aa, j, theta, knots, lambdas_new) {
  data_miss = aa[[2]]
  # the missing index of the variable Xj
  xj_miss_ind = c(1:nrow(data_miss))[is.na(data_miss[,j])]
  if(length(xj_miss_ind)==0){
    return(list(true=NULL,xj.pred=NULL))
  }
  xj.pred=matrix(0,length(xj_miss_ind),nrow(theta))
  for(i in 1:length(xj_miss_ind)){
    for(l in 1:nrow(theta)){
      xj.pred[i,l]=mu(data_miss[xj_miss_ind[i],], j, theta[l,], knots, lambdas_new[l,])
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

errors_new <- function(aa, theta, knots, lambdas_new) {
  error1 = NULL
  for (j in 1:ncol(aa[[2]])) {
    pred_true = pred_new(aa,j,theta,knots, lambdas_new)
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





#Comparison of the marginal density
#ks.test
for(j in 1:ncol(aa[[2]])){
  xpdf<-function(x){
    sum=0
    for(i in 1:nrow(b[[1]])){
      sum=sum+single_den(x,b[[1]][i,],knots[j,],lambdas[j])
    }
    return(sum/nrow(b[[1]]))
  }
  xpdff<-function(x){sapply(x,xpdf)}
  xcdf<-function(x){
    sapply(x,function(y){integrate(xpdff,-Inf,y)$value})
  }
  print(ks.test(aa[[1]][,j],"xcdf")$p)
}

xpdf<-function(x){
  sum=0
  for(i in 1:nrow(old_theta)){
    sum=sum+single_den(x,old_theta[i,],knots[j,],lambdas[j])
  }
  return(sum/nrow(old_theta))
}
xpdff<-function(x){sapply(x,xpdf)}
xcdf<-function(x){
  sapply(x,function(y){integrate(xpdff,-Inf,y)$value})
}
print(ks.test(aa[[1]][,j],"xcdf")$p)
}

#NEW WITH LAMBDA
for(j in 1:ncol(aa[[2]])){
  xpdf<-function(x){
    sum=0
    for(i in 1:nrow(c[[1]])){
      sum=sum+single_den(x,c[[1]][i,],c[[2]][j,],c[[3]][i,j])
    }
    return(sum/nrow(c[[1]]))
  }
  xpdff<-function(x){sapply(x,xpdf)}
  xcdf<-function(x){
    sapply(x,function(y){integrate(xpdff,-Inf,y)$value})
  }
  print(ks.test(aa[[1]][,j],"xcdf")$p)
}

#plot
par(mfrow=c(2,2))
for(i in 1:ncol(aa[[2]])){
  s1=summaries(aa[[1]], i, b[[1]], knots, lambdas)
  s2=summaries(aa[[1]],i,old_theta,knots,lambdas)
  x.hmax=1.2*max(hist(aa[[1]][,i],breaks=15,plot=F)$density)
  hist(aa[[1]][,i],freq=F,breaks=15,ylim=c(0,x.hmax))
  lines(s1$x.grid,s1$fx.hat.gibbs,col="blue",lwd=2)
  lines(s1$x.grid,s1$fx.hat.gibbs.l,col="blue",lty=2)
  lines(s1$x.grid,s1$fx.hat.gibbs.u,col="blue",lty=2)
  lines(s2$x.grid,s2$fx.hat.gibbs,col="red",lwd=2)
  lines(s2$x.grid,s2$fx.hat.gibbs.l,col="red",lty=2)
  lines(s2$x.grid,s2$fx.hat.gibbs.u,col="red",lty=2)
}

#plot of b and c
par(mfrow=c(2,2))
for(i in 1:ncol(aa[[2]])){
  s1=summaries(aa[[1]], i, b[[1]], knots, lambdas)
  s2=summaries_new(aa[[1]],i,c[[1]],knots_new,c[[3]])
  x.hmax=1.2*max(hist(aa[[1]][,i],breaks=15,plot=F)$density)
  hist(aa[[1]][,i],freq=F,breaks=15,ylim=c(0,x.hmax))
  lines(s1$x.grid,s1$fx.hat.gibbs,col="blue",lwd=2)
  lines(s1$x.grid,s1$fx.hat.gibbs.l,col="blue",lty=2)
  lines(s1$x.grid,s1$fx.hat.gibbs.u,col="blue",lty=2)
  lines(s2$x.grid,s2$fx.hat.gibbs,col="red",lwd=2)
  lines(s2$x.grid,s2$fx.hat.gibbs.l,col="red",lty=2)
  lines(s2$x.grid,s2$fx.hat.gibbs.u,col="red",lty=2)
}



#Comparison of the prediction of missing values
#Sum of squares
errors(aa,b[[1]],knots,lambdas)
errors(aa,old_theta,knots,lambdas)
errors_new(aa,c[[1]],knots_new,c[[3]])

#Plot
#NEW
par(mfrow=c(2,2))
for(j in 1:ncol(aa[[2]])){
  pred_true = pred(aa,j,b[[1]],knots, lambdas)
  y.pred = pred_true$xj.pred
  if(length(y.pred)!=0){
    y.pred.gibbs.l=apply(y.pred,1,quantile,prob=0.025)
    y.pred.gibbs=apply(y.pred,1,mean)
    y.pred.gibbs.u=apply(y.pred,1,quantile,prob=0.975)
    true = pred_true$true
    plot(1:length(true),sort(true),"l",col="green")
    lines(1:length(true),y.pred.gibbs[order(true)],col="blue")
#    lines(1:length(true),y.pred.gibbs.l[order(true)],col="blue",lty=2)
#    lines(1:length(true),y.pred.gibbs.u[order(true)],col="blue",lty=2)
  }
}

#NEW_LAMBDA
par(mfrow=c(2,2))
for(j in 1:ncol(aa[[2]])){
  pred_true = pred_new(aa,j,c[[1]],knots_new, c[[3]])
  y.pred = pred_true$xj.pred
  if(length(y.pred)!=0){
    y.pred.gibbs.l=apply(y.pred,1,quantile,prob=0.025)
    y.pred.gibbs=apply(y.pred,1,mean)
    y.pred.gibbs.u=apply(y.pred,1,quantile,prob=0.975)
    true = pred_true$true
    plot(1:length(true),sort(true),"l",col="green")
    lines(1:length(true),y.pred.gibbs[order(true)],col="red")
    #lines(1:length(true),y.pred.gibbs.l[order(true)],col="red",lty=2)
    #lines(1:length(true),y.pred.gibbs.u[order(true)],col="red",lty=2)
  }
}


#OLD
for(j in 1:ncol(aa[[2]])){
  pred_true = pred(aa,j,old_theta,knots, lambdas)
  y.pred = pred_true$xj.pred
  y.pred.gibbs.l=apply(y.pred,1,quantile,prob=0.025)
  y.pred.gibbs=apply(y.pred,1,mean)
  y.pred.gibbs.u=apply(y.pred,1,quantile,prob=0.975)
  true = pred_true$true
  plot(1:length(true),sort(true),"l",col="green")
  lines(1:length(true),y.pred.gibbs[order(true)],col="red")
  lines(1:length(true),y.pred.gibbs.l[order(true)],col="red",lty=2)
  lines(1:length(true),y.pred.gibbs.u[order(true)],col="red",lty=2)
}

# r package
error3 <- NULL
for (j in 1:ncol(aa[[2]])) {
  ind <- is.na(aa[[2]][,j])
  true = aa[[1]][ind,j]
  plot(1:length(true),sort(true),"l",col="green")
  imputed <- aa[[3]][ind,j][order(true)]
  lines(1:length(true),imputed,col="red")
  error3 <- c(error3, sum((sort(true)-imputed)^2))
}
error3

#Hmisc
par(mfrow=c(2,2))
error4 <- NULL
for (j in 1:ncol(aa[[2]])) {
  ind <- is.na(aa[[2]][,j])
  true = aa[[1]][ind,j]
  plot(1:length(true),sort(true),"l",col="green")
  imputed <- aa[[4]][ind,j][order(true)]
  lines(1:length(true),imputed,col="red")
  error4 <- c(error4, sum((sort(true)-imputed)^2))
}
error4

#rf.imp
par(mfrow=c(2,2))
error5 <- NULL
for (j in 1:ncol(aa[[2]])) {
  ind <- is.na(aa[[2]][,j])
  true = aa[[1]][ind,j]
  plot(1:length(true),sort(true),"l",col="green")
  imputed <- aa[[5]][ind,j][order(true)]
  lines(1:length(true),imputed,col="red")
  error5 <- c(error5, sum((sort(true)-imputed)^2))
}
error5

#am.imp
par(mfrow=c(2,2))
error6 <- NULL
for (j in 1:ncol(aa[[2]])) {
  ind <- is.na(aa[[2]][,j])
  true = aa[[1]][ind,j]
  plot(1:length(true),sort(true),"l",col="green")
  imputed <- aa[[6]][ind,j][order(true)]
  lines(1:length(true),imputed,col="red")
  error6 <- c(error6, sum((sort(true)-imputed)^2))
}
error6

#mi.imp
par(mfrow=c(2,2))
error7 <- NULL
for (j in 1:ncol(aa[[2]])) {
  ind <- is.na(aa[[2]][,j])
  true = aa[[1]][ind,j]
  plot(1:length(true),sort(true),"l",col="green")
  imputed <- aa[[7]][ind,j][order(true)]
  lines(1:length(true),imputed,col="red")
  error7 <- c(error7, sum((sort(true)-imputed)^2))
}
error7





b=density_estimate1(aa[[1]],m,100,rep(1/m,m),100)
theta=b[[1]];x.knot=b[[2]];y.knot=b[[3]];lambda.x=b[[4]];lambda.y=b[[5]]
knots=rbind(x.knot,y.knot);lambdas=c(lambda.x,lambda.y)
x=aa[[1]][,1];
x.grid=seq(min(x),max(x),l=length(x));
mu.hat=matrix(0,length(x.grid),nrow(theta))
for(i in 1:length(x.grid)){
  for(l in 1:nrow(theta)){
    mu.hat[i,l]=mu(c(x.grid[i],NA),2, theta[l,],knots,lambdas)
  }}
mu.hat.gibbs.l=apply(mu.hat,1,quantile,prob=0.025)
mu.hat.gibbs=apply(mu.hat,1,mean)
mu.hat.gibbs.u=apply(mu.hat,1,quantile,prob=0.975)

plot(x.grid,mu.hat.gibbs,col="red","l",lwd=2)
lines(x.grid,mu.hat.gibbs.l,col=2,lty=2)
lines(x.grid,mu.hat.gibbs.u,col=2,lty=2)
lines(x.grid,mu.true(x.grid),col="green",lwd=2,lty=4)


cv_miss(aa[[2]],1)
cv_miss(aa[[2]][complete.cases(aa[[2]]), ],1)





#trace of theta and lambda
#lambda
par(mfrow=c(2,2))
for(j in 1:ncol(aa[[2]])){
  plot(1:nrow(c[[3]]),c[[3]][,j],"l")
}

par(mfrow=c(4,2))
for(j in 1:ncol(b[[1]])){
  plot(1:nrow(b[[1]]),b[[1]][,j],"l")
}

par(mfrow=c(4,2))
for(j in 1:ncol(c[[1]])){
  plot(1:nrow(c[[1]]),c[[1]][,j],"l")
}

b[[3]]
c[[3]][90:101,]
density(aa[[1]][,1])$bw
density(aa[[1]][,2])$bw