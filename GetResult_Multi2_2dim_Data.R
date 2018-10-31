library(Hmisc)
library(missForest)
library(Amelia)
library(mice)
library(mi)
library(mice)
source("cv_miss_new3.r")

mu.true=function(x){exp(x/6)-x+log(x^4+1)}
sigma.true=function(x){(x^2)*exp(-abs(x))}

GetResult3 <- function(miss_prop, n){
  x=rnorm(n,0,2)
  y=mu.true(x)+sigma.true(x)*rnorm(n)
  miss_x_ind = sample(1:n, size =round(miss_prop*n),  replace = FALSE)
  miss_y_ind = sample((1:n)[-miss_x_ind], size =round(miss_prop*n),  replace = FALSE)
  x_miss = x; y_miss = y
  x_miss[miss_x_ind] = NA; y_miss[miss_y_ind] = NA
  raw_data=cbind(x,y); miss_data=cbind(x_miss,y_miss)
  
  # Get the number of NA's in each column so as to calculate the MSE
  NA_num <- apply(is.na(miss_data), 2, sum)
  
  #
  tempData <- mice(miss_data, m=5,meth='pmm',seed=500, print=F)
  completedData <- mice::complete(tempData, 1)
  #
  x_miss <- miss_data[,1]
  y_miss <- miss_data[,2]
  hmisc <- aregImpute(~x_miss+y_miss)
  x_hmisc <- apply(hmisc$imputed[[1]],1,mean)
  y_hmisc <- apply(hmisc$imputed[[2]],1,mean)
  x_miss[is.na(x_miss)] <- x_hmisc
  y_miss[is.na(y_miss)] <- y_hmisc
  #
  #rf.imp <- missForest(miss_data)$ximp
  # 
  #library(imputeR)
  #i.imp <- impute(miss_data1, lmFun="lassoR")$imp
  #
  m = 5 #m is the number of imputations, default is  5
  a.out <- amelia(miss_data, m=m) 
  a.imp <- a.out$imputations[[1]]
  for(i in 2:m) {
    a.imp <- a.imp+a.out$imputations[[i]]
  }
  a.imp <- a.imp/m
  #
  mi.out <- mi(missing_data.frame(miss_data), seed = 500)
  mi.imp <- complete(mi.out)
  nchain <- length(mi.imp)
  p = dim(raw_data)[2]
  m.imp <- mi.imp[[1]][,1:p]
  for(i in 2:nchain) 
    m.imp = m.imp+mi.imp[[i]][,1:p]
  m.imp <- as.matrix(m.imp/nchain)
  #### store raw data, missing data and imputated data from different imputation methods
  aa = list(raw = raw_data, miss = miss_data, mice.imp=as.matrix(completedData),
            hmisc.imp = cbind(x_miss,y_miss),
            am.imp = a.imp, mi.imp = m.imp)
  
  #b=density_estimate(aa[[2]],m,100,rep(1/m,m),100,0)
  
  #########NEW METHOD
  #theta=b[[1]]; knots=b[[2]]; lambdas=b[[3]];
  
  ######### NEW METHOD WITH LAMBDA
  # select m using 5 folder cv for new method
  m_5 = cv_miss_new3(aa[[2]], 5, cl)
  com_data=aa[[2]][complete.cases(aa[[2]]), ]
  p=ncol(aa[[2]])
  init2=rep(1,p)
  init1_5 = rep(1/m_5,m_5);
  for (i in 1:p) {
    y <- com_data[,i]
    init2[i] = (density(y)$bw)^2
  }
  
  d=density_estimate_new(aa[[2]],m_5,N_burn=100,init1_5,init2,N_iter=200,N_burn_out=100)
  theta_new_5=d[[1]]; knots_new_5=d[[2]]; lambdas_new_5=d[[3]]
  
  ########### OLD METHOD
  # select m using 5 folder cv for old method
  m_5_old = cv_miss_new3(com_data, 5, cl)
  init1_5_old = rep(1/m_5_old,m_5_old);
  o=density_estimate_new(com_data, m_5_old,N_burn=100,init1_5_old,init2,N_iter=200,N_burn_out=100)
  
  #Comparison of the prediction of missing values
  #ks.test
  ks_test_old <- NULL;
  ks_test_new <- NULL; 
  for(j in 1:ncol(aa[[2]])){
    xpdf<-function(x){
      sum=0
      for(i in 1:nrow(d[[1]])){
        sum=sum+single_den(x,d[[1]][i,],d[[2]][j,],d[[3]][i,j])
      }
      return(sum/nrow(d[[1]]))
    }
    xpdff<-function(x){sapply(x,xpdf)}
    xcdf<-function(x){
      sapply(x,function(y){integrate(xpdff,-Inf,y)$value})
    }
    ks_test_new <- c(ks_test_new, ks.test(aa[[1]][,j],"xcdf")$p)
    
    xpdf<-function(x){
      sum=0
      for(i in 1:nrow(o[[1]])){
        sum=sum+single_den(x,o[[1]][i,],o[[2]][j,],o[[3]][i,j])
      }
      return(sum/nrow(o[[1]]))
    }
    xpdff<-function(x){sapply(x,xpdf)}
    xcdf<-function(x){
      sapply(x,function(y){integrate(xpdff,-Inf,y)$value})
    }
    ks_test_old <- c(ks_test_old, ks.test(aa[[1]][,j],"xcdf")$p)
  }
  
  
  #Sum of squares
  # New method with lambda
  error_new_5=errors_new(aa,d[[1]],knots_new_5,d[[3]])/NA_num
  print("error_new_5")
  #print(error_new_5[iter,])
  # Old method with lambda
  error_o111=errors_new(aa,o[[1]],o[[2]],o[[3]])/NA_num
  print("error_old_5")
  #print(error_o111[iter,])
  error <- c(error_new_5, error_o111)
  for (k in 3:6) {
    error3=NULL
    for (j in 1:ncol(aa[[2]])) {
      ind <- is.na(aa[[2]][,j])
      true = aa[[1]][ind,j]
      imputed <- aa[[k]][ind,j][order(true)]
      error3 <- c(error3, sum((sort(true)-imputed)^2)/NA_num[j])
    }
    error <- c(error, error3) # length is 2*6, since each method (total 6 methods) has 2 dimensional error
  }
  err_x=error[seq(1,12,2)]
  err_y=error[seq(2,12,2)]
  ks_x <- c(ks_test_new[1], ks_test_old[1], rep(NA, 4))
  ks_y <- c(ks_test_new[2], ks_test_old[2], rep(NA, 4))

  m <- c(m_5,  m_5_old, rep(NA, 4))
  return(list(err_x=err_x, err_y=err_y, 
              ks_x=ks_x, ks_y=ks_y, m=m, n=dim(data)[1]))
}

### not use parallel
N=30;
miss_prop =0.1
for (miss_prop in c(0.1,0.2,0.4)) {
  n = 400
  err_x <- err_y <- ks_x <- ks_y <- m <- vector()
  for(i in 1:N) {
    #miss_prop=0.4
    res <- GetResult3(miss_prop, n)
    #print(res)
    err_x <- rbind(err_x, res$err_x)
    #print(err_x)
    err_y <- rbind(err_y, res$err_y)
    
    ks_x <- rbind(ks_x, res$ks_x)
    ks_y <- rbind(ks_y, res$ks_y)
    
    m <- rbind(m, res$m)
    print(paste0("Replicate: ",i))
    DATA <- rbind(err_x,err_y,ks_x,ks_y,m)
    row.names(DATA) <- rep(c("error.x","error.y",
                             "ks.pval.x","ks.pval.y","m"),rep(dim(DATA)[1]/5,5))
    colnames(DATA) <- c("New","Old","Mice","hmisc","am","mi")
    write.csv(DATA, file=paste0("\\\\wolftech.ad.ncsu.edu/cos/stat/Redirect/lwu9/Documents/740/FinalProject/newsimresults/newtwodim_n",n,"miss",miss_prop,".csv"))
  }
  
}

# #t1 <- Sys.time()
# registerDoSNOW(cl)
# errors <- foreach(i = 1:m_num,
#                   .combine = 'cbind',
#                   .export=c("density_estimate_new", 
#                             "conditional_sampler", 
#                             "posterior_Gibbs_sampler", 
#                             "rdirichlet",
#                             "posterior_Gibbs_sampler.lambda_sq")) %:% 
#   foreach(j = 1:nfolder, .combine = c) %dopar% {
#     return(squared_error(group_index[[j]], m1[i]))
#     
#   }
# #t2 <- Sys.time()
# #print(difftime(t2,t1))
# 
# N=16;n=100
# miss=c(0.1,0.2,0.4) 
# miss=0.4
# for( miss_prop in miss) {
#   err_x <- err_y <- ks_x <- ks_y <- vector()
#   r <- foreach(icount(N)) %dopar% {
#     library(mice)
#     
#     #miss_prop=0.1
#     res <- GetResult(n, miss_prop)
#     #print(res)
#     err_x <- res$err_x
#     #print(err_x)
#     err_y <- res$err_y
#     ks_x <- res$ks_x
#     ks_y <- res$ks_y
#     m <- res$m
#     rbind(err_x,err_y,ks_x,ks_y,m)
#   }
#   
#   DATA <- matrix(0, N*5, 3)
#   for (i in 1:N) {
#     for (j in 1:5)
#       DATA[(j-1)*N+i,] <- r[[i]][j,]
#     
#   }
#   row.names(DATA) <- rep(c("error.x","error.y", "ks.pval.x","ks.pval.y","m"),rep(N,5))
#   colnames(DATA) <- c("New","Old","Mice")
#   write.csv(DATA, file=paste0("S:\\Documents\\740\\FinalProject\\simreults\\twodim2_n_new",n,"miss",miss_prop,".csv"))
#   
# }
# 
# 
