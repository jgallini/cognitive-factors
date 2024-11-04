
#Forward filtering backward sampling function.
source("ffbsUneq2.R")

#Data needs to group by id, then sorted by time

BayesKalm.Uneq <- function(y.long, X.long, id, time, Burn = 500, Its = 1500, Beta.Initial = 0, sigma2.beta = 10, u0 = 0, P0 = 10, a0 = 0.01, b0 = 0.01, c0 = 0.01, d0 = 0.01, silence = FALSE){
  p <- ncol(X.long)
  #hard coding this for now, used to be: ncol(X.long)
  n <- length(unique(id))
  n.totobs <- length(y.long)
  T <- max(table(id))
  
  
  ###CONVERT FROM LONG TO SHORT DATA
  #Configure the outcome
  yl <- lapply(unique(id), function(x){
    ys <- y.long[id == x]
    c(ys, rep(NA, T - length(ys)))
  })
  y <- do.call("rbind", yl)
  isnay <- y ^ is.na(y)
  #Configure time   
  tl <- lapply(unique(id), function(x){
    times <- time[id == x]
    c(times, rep(NA, T - length(times)))
  })
  timeMat <- do.call("rbind", tl)
  timeDiff <- apply(timeMat, 1, diff)
  timeDiff[is.na(timeDiff)] <- 0
  
  NON.MISSING <- !is.na(y)
  cat("check")
  #Configure X
  
 # xl <- lapply(unique(id), function(x){
#    xs <- X.long[id == x,]
#    rbind(xs, matrix(0, nrow = T-nrow(xs), ncol = 3))
#  })
  
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
  if(length(Beta.Initial)  == p) Beta.Initial <- Beta.Initial else Beta.Initial <- rep(Beta.Initial, p)
  B.star <- Beta.Initial
  
  # initialize starting points for the variance parameters
  sigma2.eta.star <- 1
  sigma2.eps.star <- 1
  
  # Create empty objects to fill in KF and KS
  Beta.Track <- matrix(NA, p, Its)
  sigma2.eps.Track <- numeric(Its)
  sigma2.eta.Track <- numeric(Its)
  y.star.Track <- mu.Track <- array(NA, dim = c(n, T, Its))
  
  ######################################################################################
  nv.track <- numeric(Its)
  nv.length <- numeric(Its)
  
  #Needed for beta posterior
  sig2beta_XtX <- sigma2.beta * crossprod(X.long)
  ex <- eigen(sig2beta_XtX, symm = TRUE)
  
  
  if(!silence) pb <- txtProgressBar(min = 0, max = Its, style = 3)
  
  for(j in 1:Its){
    
    ###### mu post
    y.star <- y - apply(X,3, function(x) x %*% B.star)
    y.star.Track[,,j] <- y.star    
    
    mu.star <- suppressWarnings(ffbs.uneq2(y = y.star, V = sigma2.eps.star, W = sigma2.eta.star, m0 = u0, C0 = P0, timeDiff = timeDiff)$x)
    mu.Track[,,j] <- mu.star
    
    
    ###### sig 2 eta Post
    
    nu <- sum(apply(mu.star * isnay, 1, diff)^2/timeDiff, na.rm = TRUE)
    sigma2.eta.star <- sigma2.eta.Track[j] <- 1/rgamma(1, ((n.totobs-n)/ 2 + a0), (b0 + nu/ 2))
    
    
    ##### sig 2 eps Post
    nv <- sum((y.star-mu.star)^2, na.rm = TRUE)
    sigma2.eps.star <- sigma2.eps.Track[j] <- 1/rgamma(1, (n.totobs/2 +c0), d0+nv/2)
    
    
    # #### Beta Post
    
    
    v.star <- (y-mu.star)
    B.sum <- rowSums(sapply(1:T, function(x) colSums(v.star[,x] * X[,,x], na.rm = TRUE)), na.rm = TRUE)
    B.Big <- sigma2.beta * B.sum + sigma2.eta.star * Beta.Initial
    Sigma.Inv <- tcrossprod(ex$vectors/(ex$values + sigma2.eps.star)[col(ex$vectors)], ex$vectors)
    B.star <- Beta.Track[,j] <- mvtnorm::rmvnorm(1, mean = crossprod(Sigma.Inv, B.Big), sigma = Sigma.Inv * sigma2.eps.star * sigma2.beta)[1,]
    
    if(!silence) setTxtProgressBar(pb, j)
  }
  cat("made it here")
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
  assign("bkout",bkout,envir=.GlobalEnv)
  cat("made it even further")
}

































