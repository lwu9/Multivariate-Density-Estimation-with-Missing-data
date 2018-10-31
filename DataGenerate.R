 
library(mice)
data <- airquality
dataT <- data[complete.cases(data), 1:4]
data <- dataT
data[4:10,3] <- rep(NA,7)
data[1:5,4] <- NA
tempData <- mice(data,m=5,maxit=50,meth='pmm',seed=500)
summary(tempData)
completedData <- complete(tempData,1)
aa <- list(true=as.matrix(dataT), miss=as.matrix(data), imputed=as.matrix(completedData))

## airquality
data(airquality)
head(airquality)
data <- airquality
dataT <- data[complete.cases(data), 1:4]
data <- dataT
x = data[,1]; y=data[,2]; z=data[,3]; q=data[,4]
n = nrow(data)
miss_x_prop = 0.1; miss_y_prop= 0.1; miss_z_prop=0.1; miss_q_prop=0.1
miss_x_ind = sample(1:n, size =round(miss_x_prop*n),  replace = FALSE)
miss_y_ind = sample(1:n, size =round(miss_y_prop*n),  replace = FALSE)
miss_z_ind = sample(1:n, size =round(miss_z_prop*n),  replace = FALSE)
miss_q_ind = sample(1:n, size =round(miss_z_prop*n),  replace = FALSE)
x_miss = x; y_miss = y; z_miss=z; q_miss=q
x_miss[miss_x_ind] = NA; y_miss[miss_y_ind] = NA; z_miss[miss_z_ind] = NA; q_miss[miss_q_ind] = NA
raw_data=cbind(x,y,z,q); miss_data=cbind(x_miss,y_miss,z_miss,q_miss)
raw_data1=miss_data1=NULL
for(i in 1:n){
  if(sum(is.na(miss_data[i,]))!=4){
    raw_data1=rbind(raw_data1,raw_data[i,])
    miss_data1=rbind(miss_data1,miss_data[i,])
  }
}
tempData <- mice(miss_data1, m=5, maxit=50,meth='pmm',seed=500)
summary(tempData)
completedData <- complete(tempData,1)
aa = list(raw_data1, miss_data1, imputed=as.matrix(completedData))


b=density_estimate(aa[[2]],m,100,rep(1/m,m),100)

# Simulate data set with \code{mvrnorm} from package \code{\pkg{MASS}}.
require(MASS)
sigma <- matrix(data = c(1, 1.2, 2.2, 1.2, 1, 1.2, 2.2, 1.2, 1), nrow = 3)
sigma <- sigma%*%t(sigma)
complete.data <- mvrnorm(n = 100, mu = c(5, 5, 5), Sigma = sigma)
# Perform quick amputation
result1 <- ampute(data = complete.data)
# Change default matrices as desired
patterns <- result1$patterns
patterns[1:3, 2] <- 0
odds <- result1$odds
odds[2,3:4] <- c(2, 4)
odds[3,] <- c(3, 1, NA, NA)
# Rerun amputation
result3 <- ampute(data = complete.data, patterns = patterns, freq = 
                    c(0.3, 0.3, 0.4), cont = FALSE, odds = odds)
# Run an amputation procedure with continuous probabilities
result4 <- ampute(data = complete.data, type = c("RIGHT", "TAIL", "LEFT"))

# }
tempData <- mice(result1$amp,m=5,maxit=50,meth='pmm',seed=500)
completedData <- complete(tempData,1)
aa <- list(true=as.matrix(complete.data), miss=as.matrix(result1$amp), imputed=as.matrix(completedData))

## data in class
#Generate a raw data set and a missing data set
set.seed(740)
#set true mean and sd functions
mu.true=function(x){exp(x/6)-x+log(x^4+1)}
sigma.true=function(x){(x^2)*exp(-abs(x))}
#par(mfrow=c(2,2))
#curve(mu.true,-5,5); curve(sigma.true,-5,5)

#generate data:
#set.seed(740)
n=100; x=rnorm(n,0,2)
y=mu.true(x)+sigma.true(x)*rnorm(n)
miss_x_prop = 0.4; miss_y_prop= 0.4
miss_x_ind = sample(1:n, size =round(miss_x_prop*n),  replace = FALSE)
miss_y_ind = sample((1:n)[-miss_x_ind], size =round(miss_y_prop*n),  replace = FALSE)
x_miss = x; y_miss = y
x_miss[miss_x_ind] = NA; y_miss[miss_y_ind] = NA
tempData <- mice(cbind(x_miss, y_miss),m=5,maxit=50,meth='pmm',seed=500)
completedData <- complete(tempData,1)
aa = list(cbind(x,y), cbind(x_miss, y_miss), imputed = as.matrix(completedData))

