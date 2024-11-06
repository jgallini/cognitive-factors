
# This function runs the forward Kalman filter and the backward smoother
# to estimate alpha values within the Gibbs sampler in BayesKalmJoint_gh.R

# Inputs:
# y (a matrix of outcomes)
# V (a diagonal matrix of initial variance estimates for Sigma epsilon)
# W (a diagonal matrix of initial variance estimates for Sigma eta)
# m0 (scalar initial estimate for alpha mean)
# C0 (scalar initial estimate for alpha variance)
# timeDiff (matrix of time difference between visits)

ffbs.joint = function(y,V,W,m0,C0, timeDiff){
  
  #Initializing values and creating storage matrices
  T = ncol(y); n = nrow(y)
  a = matrix(0, n, T); R = array(0, dim = c(n, n, T)) 
  m = matrix(0, n, T); C = array(0, dim = c(n, n, T)) ; 
  B = array(0, dim = c(n, n, T-1))
  H = array(0, dim = c(n, n, T-1)); mm = matrix(0, n, T); 
  CC = array(0, dim = c(n, n, T))
  x = matrix(0, n, T); llike = 0.0 
  
  #Beginning forward Kalman filter calculations
  for (t in 1:T){ 
    if(t==1){
      a[,1] = m0; R[,,1] = C0 + W
    }else{
      a[,t] = m[,t-1]; R[,,t] = C[,,t-1] + timeDiff[t-1] * W
    }
    f = a[,t]
    Q = R[,,t] + V 
    A = R[,,t] %*% solve(Q)
    Av = A %*% (y[,t]-f)
    m[,t] = a[,t]+Av 
    C[,,t] = R[,,t]- A %*% R[,,t]
    
    if(t > 1){
      B[,,t-1] = C[,,t-1] %*% solve(R[,,t])
      H[,,t-1] = C[,,t-1]- B[,,t - 1] %*% tcrossprod(R[,,t],B[,,t - 1])
    }
  }
  
  #Prepping values for backward sampler
  mm[,T] = m[,T]; CC[,,T] = C[,,T]
  x[,T] = rmvnorm(1, mean = m[,T], sigma = (C[,,T]), checkSymmetry = FALSE) 
  
  #Beginning backward sampler
  for (t in (T-1):1){ 
    x[,t] = rmvnorm(1,mean = m[,t]+B[,,t]%*%(x[,t+1]-a[,t+1]), 
                    sigma = H[,,t], checkSymmetry = FALSE)
  } 
  #storing estimates of interest
  return(list(x=x,m=m,C=C,mm=mm,CC=CC,llike=llike))
}
