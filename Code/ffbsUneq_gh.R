
# This function runs the forward Kalman filter and the backward sampler to 
# calculate estimates for the mean and variance of the alpha parameters.
# This particular ffbs code is used for estimating the alpha
# values in a univariate state-space model (one outcome at a time).
# In context of the multivariate state-space model, this code is 
# used to estimate the initial values for alpha to begin the Gibbs
# sampler that estimates the multivariate model.

# Inputs:
# y (matrix of outcome data)
# V (scalar variance estimate of sigma epsilon)
# W (scalar variance estimate of sigma eta)
# m0 (scalar initial alpha estimate)
# C0 (scalar initial variance of alpha estimate)
# timeDiff (matrix of differences in time between each visit)

ffbs.uneq2 = function(y,V,W,m0,C0, timeDiff){ 
  
  #Initializing values and creating storage matrices
  T = ncol(y); n = nrow(y)
  #Predicted alpha storage matrix
  a = matrix(0, n, T)
  #Predicted alpha covariance storage matrix
  R = matrix(0, n, T)
  #Updated alpha estimates storage matrix
  m = matrix(0, n, T)
  #Updated alpha covariance storage matrix
  C = matrix(0, n, T)
  #storage for derived matrix needed in backward smoother
  B = matrix(0, n, T-1)
  #storage for covariance matrix for backward smoother
  H = matrix(0, n, T) #changing this to be dimension T instead of T-1
  #storage for final alpha estimates, this will be main output of this function
  x = matrix(0, n, T)
  #additional values of possible interest
  mm = matrix(0, n, T) 
  CC = matrix(0, n, T)
  llike = 0.0 
  is.na.y <- is.na(y)
  #initiating time counting variables
  t <- 0
  t <- t+1
  
  #beginning forward Kalman filter calculations
  for (t in 1:T){ 
    if(t==1){
      #initiate alpha estimate and covariance of alpha estimate
      a[,1] = m0; R[,1] = C0 + W
    }else{
      #take most recent updated estimates for alpha and covariance of alpha
      a[,t] = m[,t-1]; R[,t] = C[,t-1] + timeDiff[t-1, ] * W
    }
    #calculating updated alpha estimate and covariance
    #estimate using Kalman filter equations
    f = a[,t]
    Q = R[,t] + V 
    #A is Kalman gain
    A = R[,t]/Q
    Av = ifelse(is.na(y[,t]), 0, A * (y[,t]-f))
    #updated alpha estimate
    m[,t] = a[,t]+Av 
    #updated alpha covariance estimate
    C[,t] = R[,t]-Q*A**2 
    #additional values of interest for smoother
    B[,t-1] = C[,t-1]/R[,t] 
    llike = llike + sum(dnorm(y[,t],f,sqrt(Q),log=TRUE), na.rm = TRUE) 
  }
  
  #prepping values for backward sampler
  mm[,T] = m[,T]; CC[,T] = C[,T]
  C <- C * !is.na.y
  B <- B * !is.na.y[,-1]
  #randomly drawing alpha estimates to initiate smoother from normal
  #distribution defined by updated estimates and updated covariance
  x[,T] = rnorm(n,m[,T],sqrt(C[,T])) 
  #initiating H[,t] to be equal to updated variance estimate at time T
  H[,T] = C[,T]
  #starting time at t-1
  t <- T-1
  
  #beginning smoothing backward sampler calculations
  for (t in (T-1):1){ 
    #smoothed predicted mean alpha estimate, can output if interested in this
    mm[,t] = m[,t] + C[,t]/R[,t+1]*(mm[,t+1]-a[,t+1])
    #smoothed predicted alpha covariance estimate, 
    #can output if interested in this
    CC[,t] = C[,t] - (C[,t]^2)/(R[,t+1]^2)*(R[,t+1]-CC[,t+1])
    #calculating smoothed variance
    H[,t] = C[,t]-(H[,t+1]-R[,t+1])*B[,t]**2
    #random draw from smoothed alpha posterior distribution, 
    #main estimates of interest
    x[,t] = rnorm(n,m[,t]+B[,t]*(x[,t+1]-a[,t+1]),
                  sqrt(H[,t] + C[,t] * (H[,t] == 0)))
    #just in case H[,t] is 0, shouldn't be in most cases
  } 
  #returning estimates of interest, x is the main alpha estimates of interest
  return(list(x=x,m=m,C=C,mm=mm,CC=CC,llike=llike)) 
}

