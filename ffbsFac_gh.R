
#This program runs the Kalman filter and backward smoother to estimate the 
#alpha values in the model
#It runs for each person one by one

# Inputs:
# y: a matrix of the outcome variables for all time points for one person
# G: a matrix that represents the most recent G estimate in the sampler
# V: a matrix that represents the most recent sigma epsilon estimate in
#the sampler
# W: a matrix that represents the most recent sigma eta estimate in the 
#sampler
# m0: a scalar that represents the initial starting value for each alpha
#value, this parameter is pretty arbitrary but recommended to use 0
# C0: a QxQ diagonal matrix that represents the initial value for the 
#covariance matrix of the alphas, also relatively arbitrary recommended to
# use a diagonal matrix of 10's
# timeDiff: a vector of the differences in time between each observation
#for one person, defaults to a vector of 1's if not overwritten

#Function outputs: alpha estimates and covariance matrices

ffbs.fac = function(y,G,V,W,m0,C0, timeDiff = rep(1, ncol(y)-1)){ #should't this be ncol(y)-1????
  
  #Variable initialization
  T = ncol(y); n = nrow(W)
  #alpha storage matrix
  a = matrix(0, n, T)
  #Predicted alpha covariance storage matrix
  R = array(0, dim = c(n, n, T)) 
  #Updated alpha estimates storage matrix
  m = matrix(0, n, T)
  #Updated alpha covariance storage matrix
  C = array(0, dim = c(n, n, T))
  #storage for derived matrix needed in backward smoother
  B = array(0, dim = c(n, n, T-1))
  #storage for covariance matrix for backward smoother
  H = array(0, dim = c(n, n, T)) 
  #storage for final alpha estimates, this will be main output of this function
  x = matrix(0, n, T)
  
  #Begin forward Kalman filter
  t <- 1
  for (t in 1:T){ 
    if(t==1){
      #initiate alpha estimate and covariance of alpha estimate
      a[,1] = m0; R[,,1] = C0 + W 
    }else{
      #take most recent smoothed estimates for alpha and covariance of alpha
      a[,t] = m[,t-1]; R[,,t] = C[,,t-1] + timeDiff[t-1] * W
    }
    #calculating updated alpha estimate and covariance
    #estimate using Kalman filter equations
    f = a[,t]
    RG = tcrossprod(R[,,t], G)
    Q = G %*% RG + V 
    A = RG %*% solve(Q) #A is Kalman gain
    Av = A %*% (y[,t]-G %*% f)
    #updated alpha estimate
    m[,t] = a[,t]+Av
    #updated alpha covariance estimate
    C[,,t] = R[,,t]- A %*% G %*% R[,,t]
    #calculating matrix needed in backward smoother 
    if(t > 1){
      B[,,t-1] = C[,,t-1] %*% solve(R[,,t])
    }
    
  }
  
  #Drawing randomly from normal posterior distribution defined by 
  #updated alpha estimate and updated covariance of alpha
  #Setting this estimate as alpha of last time point T
  #heading into the backward smoother
  x[,T] = mvtnorm::rmvnorm(1, mean = m[,T], sigma = (C[,,T]), 
                           checkSymmetry = FALSE) 
  
  #initializing H as updated alpha variance/covariance estimate for 
  #last time point T to be used in backward smoother
  H[,,T] = C[,,T] 
  
  #Beginning backward Kalman smoother/sampler
  for (t in (T-1):1){ 
    #calculating smoothed variance/covariance using backward smoothing
    #equations
    H[,,t] = C[,,t]-(B[,,t]%*%(H[,,t+1]-R[,,t+1])%*%t(B[,,t]))
    
    #randomly drawing final alpha estimate from normal distribution defined
    #by smoothed mean and smoothed variance
    x[,t] = mvtnorm::rmvnorm(1,mean = m[,t]+B[,,t]%*%(x[,t+1]-a[,t+1]), 
                             sigma = H[,,t], checkSymmetry = FALSE)
  } 
  
 #returning alpha estimates (x) and covariance estimates (m,C) 
  return(list(x=x,m=m,C=C))
}
