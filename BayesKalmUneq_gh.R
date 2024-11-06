
# This function calculates initial estimates for all of the parameters 
# needed in the multivariate model. These estimates are then used 
# to initiate the main Gibbs sampler in BayesKalmJoint_gh.R.

# Inputs
# y.long (a long matrix of outcome data, in this case one outcome)
# X.long (a long matrix of all predictors)
# id (the ID variable to be used for each participant)
# time (the time variable to be used)
# Burn (scalar number of burn-in iterations, default 500)
# Its (scaler number of iterations, default=1500)
# Beta.Initial (initial scalar value for beta, default=0)
# sigma2.beta (initial scalar value for variance of beta, default=10) 
# u0 (initial scalar value for alpha, default=0)
# P0 (initial scalar value for variance of alpha, default=10)
# a0,b0,c0,d0 (scalar prior parameters for sigma eta and sigma epsilon, 
# default = 0.01)
# silence (TRUE/FALSE value for whether or not to suppress update output, 
# default=FALSE)

#Loading Kalman filter function.
source("ffbsUneq_gh.R")

BayesKalm.Uneq <- function(y.long, X.long, id, time, Burn = 500, 
                           Its = 1500, Beta.Initial = 0, sigma2.beta = 10, 
                           u0 = 0, P0 = 10, a0 = 0.01, b0 = 0.01, c0 = 0.01, 
                           d0 = 0.01, silence = FALSE){
  
  #calculating values needed in later calculations
  p <- ncol(X.long)
  n <- length(unique(id))
  n.totobs <- length(y.long)
  T <- max(table(id))
  
  #Convert the outcome from long to short data, pad with NA's
  #where participants don't have the max number of measurements
  yl <- lapply(unique(id), function(x){
    ys <- y.long[id == x]
    c(ys, rep(NA, T - length(ys)))
  })
  y <- do.call("rbind", yl)
  #creating logical matrix indicating whether an outcome is missing or not
  isnay <- y ^ is.na(y)
  
  #Configure times in a similar fashion to outcomes above  
  tl <- lapply(unique(id), function(x){
    times <- time[id == x]
    c(times, rep(NA, T - length(times)))
  })
  timeMat <- do.call("rbind", tl)
  timeDiff <- apply(timeMat, 1, diff)
  timeDiff[is.na(timeDiff)] <- 0
  
  NON.MISSING <- !is.na(y)
  
  #printing check while simulation is running
  cat("check")
  
  #Configure X similarly to y and time above
  xl <- lapply(unique(id), function(x){
    xs <- data.matrix(X.long[id == x,])
    zeros = matrix(0, nrow = T-nrow(xs), ncol = ncol(xs))
    colnames(zeros) = colnames(xs)
    rbind(xs, zeros)
  })
  X <- array(NA, dim = c(n, p, T))
  for(i in 1:length(xl))
    X[i,,] <- t(xl[[i]])
  
  # initialize beta
  if(length(Beta.Initial)  == p) Beta.Initial <- Beta.Initial 
  else Beta.Initial <- rep(Beta.Initial, p)
  B.star <- Beta.Initial
  
  # initialize starting points for the variance parameters
  sigma2.eta.star <- 1
  sigma2.eps.star <- 1
  
  # Create empty objects to fill in during estimation
  Beta.Track <- matrix(NA, p, Its)
  sigma2.eps.Track <- numeric(Its)
  sigma2.eta.Track <- numeric(Its)
  y.star.Track <- mu.Track <- array(NA, dim = c(n, T, Its))
  
  nv.track <- numeric(Its)
  nv.length <- numeric(Its)
  
  #Calculations needed for beta posterior
  sig2beta_XtX <- sigma2.beta * crossprod(X.long)
  ex <- eigen(sig2beta_XtX, symm = TRUE)
  
  #printing update to see while code is running
  if(!silence) pb <- txtProgressBar(min = 0, max = Its, style = 3)
  
  #beginning Gibbs sampler
  for(j in 1:Its){
    
    #Alpha estimation (called mu in some cases)
    #prepping y data by subtracting linear effects
    y.star <- y - apply(X,3, function(x) x %*% B.star)
    y.star.Track[,,j] <- y.star    
    #getting alpha estimates from Kalman filter
    mu.star <- suppressWarnings(ffbs.uneq2(y = y.star, V = sigma2.eps.star, 
                                           W = sigma2.eta.star, m0 = u0, 
                                           C0 = P0, timeDiff = timeDiff)$x)
    mu.Track[,,j] <- mu.star
    
    
    #Sigma Eta estimation
    nu <- sum(apply(mu.star * isnay, 1, diff)^2/timeDiff, na.rm = TRUE)
    #Drawing from the posterior of sigma eta and storing as current estimate
    sigma2.eta.star <- 
      sigma2.eta.Track[j] <- 
      1/rgamma(1, ((n.totobs-n)/ 2 + a0), (b0 + nu/ 2))
    
    #Sigma Epsilon estimation
    nv <- sum((y.star-mu.star)^2, na.rm = TRUE)
    #Drawing from posterior of sigma epsilon and storing as current estimate
    sigma2.eps.star <- 
      sigma2.eps.Track[j] <- 1/rgamma(1, (n.totobs/2 +c0), d0+nv/2)
    
    
    #Beta estimation
    v.star <- (y-mu.star)
    B.sum <- rowSums(sapply(1:T, function(x) colSums(v.star[,x] * X[,,x], 
                                                na.rm = TRUE)), na.rm = TRUE)
    B.Big <- sigma2.beta * B.sum + sigma2.eta.star * Beta.Initial
    Sigma.Inv <- tcrossprod(
      ex$vectors/(ex$values + sigma2.eps.star)[col(ex$vectors)], ex$vectors)
    #Drawing from posterior of beta and storing as current estimate
    B.star <- Beta.Track[,j] <- 
      mvtnorm::rmvnorm(1, mean = crossprod(Sigma.Inv, B.Big), 
                       sigma = Sigma.Inv * sigma2.eps.star * sigma2.beta)[1,]
    
    #printing progress estimate while code is running
    if(!silence) setTxtProgressBar(pb, j)
  }
  
  #storing list of estimates of interest
  bkout<-list(
    Beta = Beta.Track,
    sigma2.eps = sigma2.eps.Track,
    sigma2.eta = sigma2.eta.Track,
    mu = mu.Track,
    X = X,
    timeMat = timeMat,
    y = y, 
    id = unique(id)
  )
  
  #storing output globally
  assign("bkout",bkout,envir=.GlobalEnv)
}

































