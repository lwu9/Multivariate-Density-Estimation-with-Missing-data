library(Hmisc)
library(missForest)
library(Amelia)
library(mice)
library(mi)
library(mice)

it=30
error_new_5=matrix(0,it,4)
error_new_10=matrix(0,it,4)
error_hmisc=matrix(0,it,4)
error_rf=matrix(0,it,4)
error_am=matrix(0,it,4)
error_mi=matrix(0,it,4)

for(iter in 1:it){
  
  print(iter)
  
  data(airquality)
  head(airquality)
  data <- airquality
  dataT <- data[complete.cases(data), 1:4]
  data <- dataT
  x = data[,1]; y=data[,2]; z=data[,3]; q=data[,4]
  n = nrow(data)
  miss_x_prop = 0.4; miss_y_prop= 0.4; miss_z_prop=0.4; miss_q_prop=0.4
  miss_x_ind = sample(1:n, size =round(miss_x_prop*n),  replace = FALSE)
  miss_y_ind = sample(1:n, size =round(miss_y_prop*n),  replace = FALSE)
  miss_z_ind = sample(1:n, size =round(miss_z_prop*n),  replace = FALSE)
  miss_q_ind = sample(1:n, size =round(miss_z_prop*n),  replace = FALSE)
  x_miss = x; y_miss = y; z_miss=z; q_miss=q
  x_miss[miss_x_ind] = NA; y_miss[miss_y_ind] = NA; z_miss[miss_z_ind] = NA; q_miss[miss_q_ind] = NA
  raw_data=cbind(x,y,z,q); miss_data=cbind(x_miss,y_miss,z_miss,q_miss)
  raw_data1=miss_data1=NULL
  # delete the rows which are all NA's
  for(i in 1:n){
    if(sum(is.na(miss_data[i,]))!=4){
      raw_data1=rbind(raw_data1,raw_data[i,])
      miss_data1=rbind(miss_data1,miss_data[i,])
    }
  }
  #tempData <- mice(miss_data1, m=5,meth='pmm',seed=500)
  #completedData <- complete(tempData,1)
  
  #
  x_miss <- miss_data1[,1]
  y_miss <- miss_data1[,2]
  z_miss <- miss_data1[,3]
  q_miss <- miss_data1[,4]
  hmisc <- aregImpute(~x_miss+y_miss+z_miss+q_miss)
  x_hmisc <- apply(hmisc$imputed[[1]],1,mean)
  y_hmisc <- apply(hmisc$imputed[[2]],1,mean)
  z_hmisc <- apply(hmisc$imputed[[3]],1,mean)
  q_hmisc <- apply(hmisc$imputed[[4]],1,mean)
  x_miss[is.na(x_miss)] <- x_hmisc
  y_miss[is.na(y_miss)] <- y_hmisc
  z_miss[is.na(z_miss)] <- z_hmisc
  q_miss[is.na(q_miss)] <- q_hmisc
  #
  rf.imp <- missForest(miss_data1)$ximp
  # 
  #library(imputeR)
  #i.imp <- impute(miss_data1, lmFun="lassoR")$imp
  #
  m = 5 #m is the number of imputations, default is  5
  a.out <- amelia(miss_data1, m=m) 
  a.imp <- a.out$imputations[[1]]
  for(i in 2:m) {
    a.imp <- a.imp+a.out$imputations[[i]]
  }
  a.imp <- a.imp/m
  #
  mi.out <- mi(missing_data.frame(miss_data1), seed = 500)
  mi.imp <- complete(mi.out)
  nchain <- length(mi.imp)
  p = dim(raw_data1)[2]
  m.imp <- mi.imp[[1]][,1:p]
  for(i in 2:nchain) 
    m.imp = m.imp+mi.imp[[i]][,1:p]
  m.imp <- as.matrix(m.imp/nchain)
  #### compare different methods
  aa = list(raw = raw_data1, miss = miss_data1, mice.imp=raw_data1,hmisc.imp = cbind(x_miss,y_miss,z_miss,q_miss),
            rf.imp = rf.imp, am.imp = a.imp, mi.imp = m.imp)
  
  
  
  m_5=cv_miss_new(aa[[2]],5)
  m_10=cv_miss_new(aa[[2]],10)
  print("m_5")
  print(m_5)
  print("m_10")
  print(m_10)
  
  #b=density_estimate(aa[[2]],m,100,rep(1/m,m),100,0)
  
  #########NEW METHOD
  #theta=b[[1]]; knots=b[[2]]; lambdas=b[[3]];
  
  #########NEW METHOD WITH LAMBDA
  com_data=aa[[2]][complete.cases(aa[[2]]), ]
  p=ncol(aa[[2]])
  init2=rep(1,p)
  init1_5 = rep(1/m_5,m_5);
  init1_10 = rep(1/m_10,m_10);
  for (i in 1:p) {
    y <- com_data[,i]
    init2[i] = density(y)$bw^2
  }
  
  c=density_estimate_new(aa[[2]],m_5,N_burn=100,init1_5,init2,N_iter=100,N_burn_out=0)
  theta_new=c[[1]]; knots_new=c[[2]]; lambdas_new=c[[3]]
  
  d=density_estimate_new(aa[[2]],m_10,N_burn=100,init1_10,init2,N_iter=100,N_burn_out=0)
  theta_new_10=d[[1]]; knots_new_10=d[[2]]; lambdas_new_10=d[[3]]
  
  
  
  
  #Comparison of the prediction of missing values
  #Sum of squares
  error_new_5[iter,]=errors_new(aa,c[[1]],knots_new,c[[3]])
  print("error_new_5")
  print(error_new_5[iter,])
  error_new_10[iter,]=errors_new(aa,d[[1]],knots_new_10,d[[3]])
  print("error_new_10")
  print(error_new_10[iter,])
  
  
  #hmisc
  error3=NULL
  for (j in 1:ncol(aa[[2]])) {
    ind <- is.na(aa[[2]][,j])
    true = aa[[1]][ind,j]
    imputed <- aa[[4]][ind,j][order(true)]
    error3 <- c(error3, sum((sort(true)-imputed)^2))
  }
  error_hmisc[iter,]=error3
  print("error_hmisc")
  print(error_hmisc[iter,])
  
  #rf
  error3=NULL
  for (j in 1:ncol(aa[[2]])) {
    ind <- is.na(aa[[2]][,j])
    true = aa[[1]][ind,j]
    imputed <- aa[[5]][ind,j][order(true)]
    error3 <- c(error3, sum((sort(true)-imputed)^2))
  }
  error_rf[iter,]=error3
  print("error_rf")
  print(error_rf[iter,])
  
  #am
  error3=NULL
  for (j in 1:ncol(aa[[2]])) {
    ind <- is.na(aa[[2]][,j])
    true = aa[[1]][ind,j]
    imputed <- aa[[6]][ind,j][order(true)]
    error3 <- c(error3, sum((sort(true)-imputed)^2))
  }
  error_am[iter,]=error3
  print("error_am")
  print(error_am[iter,])
  
  #mi
  error3=NULL
  for (j in 1:ncol(aa[[2]])) {
    ind <- is.na(aa[[2]][,j])
    true = aa[[1]][ind,j]
    imputed <- aa[[4]][ind,j][order(true)]
    error3 <- c(error3, sum((sort(true)-imputed)^2))
  }
  error_mi[iter,]=error3
  print("error_mi")
  print(error_mi[iter,])
  
}



apply(ks_b111,2,mean)
apply(ks_c111,2,mean)
apply(ks_o111,2,mean)

apply(error_new,2,mean)
apply(error_old,2,mean)
apply(error_hmisc,2,mean)
apply(error_rf,2,mean)
apply(error_am,2,mean)
apply(error_mi,2,mean)
