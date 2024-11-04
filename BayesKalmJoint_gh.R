
library(mvtnorm)
library(tidyverse)
source("BayesKalmUneq2_Julia.R")
source("ffbsJoint.R")

BayesKalmJoint <- function(data, outcomes, predictors, 
                                 timevar, id, initialization = "Bayes", 
                                 numits = 1000, silence = FALSE, 
                                 seed = NULL, numitInit = 1500, 
                                 burnInit = 500) {
  
  cs <- length(outcomes)
  
  if(!is.null(seed))set.seed(seed)
  
  # DataFormat ------------------------------------------------------------
  #Nest data for each id
  
  data[["time"]] <- data[[timevar]]
  data[["id"]] <- data[[id]]
  n <- length(unique(data$id))
  ndat <- data %>%
    group_by(id)%>%
    mutate(timeDiff = c(diff(time), NA)) %>%
    nest() %>%
    .$data 
  #Define matrix of output for each id
  y <- ndat %>%
    map(~select(.x, all_of(outcomes))) %>%
    map(~t(as.matrix(.x)))
  #Define matrix of predictors for each id
  X <- ndat %>%
    map(~select(.x, all_of(predictors))) %>%
    map(~(as.matrix(.x)))
  #create vector of time and timediff for each id
  time <- ndat %>%
    map("time")
  
  timeDiff <- ndat %>%
    map("timeDiff") %>%
    map(~.x[!is.na(.x)])
  
  
  
  # bayesian initialization -------------------------------------------------
  
  if(initialization == "Bayes"){
    
    BayesInit <- map(outcomes, function(outcomes){

      initmodel <- lm(as.formula(paste(outcomes, "~", paste(predictors, collapse = "+"))), data = data)
      assign("initmodel",initmodel,envir=.GlobalEnv)
      
      data$resid <- resid(initmodel)
      
      Beta.Initial <- coef(initmodel)[-1]
      
      u0 <- coef(initmodel)[1]
      
      P0 <- data %>%
        group_by(id) %>%
        filter(time == min(time)) %>%
        .[["resid"]] %>%
        var()
      
      BayesKalm.Uneq_Julia(
        y.long = data[[outcomes]], X.long = as.matrix(data[predictors]),
        id = data[[id]], time = data[[timevar]],
        Burn = burnInit, Its = numitInit, 
        Beta.Initial = Beta.Initial, sigma2.beta = 10, 
        u0 = u0, P0 = P0, 
        a0 = 0.01, b0 = 0.01, c0 = 0.01, d0 = 0.01, 
        silence = FALSE
      )
    
      #mynames <<- names(bkout)
      cat('FINISHED THISSSSSS')
      list(
        Beta.Initial <- apply(bkout$Beta[,burnInit:numitInit], 1, mean), 
        sigma2.eta <- mean(bkout$sigma2.eta[burnInit:numitInit]),
        sigma2.eps <- mean(bkout$sigma2.eps[burnInit:numitInit]),
        u0 <- u0,
        P0 <- P0,
        bkout <- list(Beta = bkout$Beta, sigma2.eta = bkout$sigma2.eta, sigma2.eps = sigma2.eps)
      )
     
      
    })
    assign("BayesInit",BayesInit,envir=.GlobalEnv)
  }
  #Julia manually setting names since they aren't taking above for some reason
BayesInit<-lapply(BayesInit,setNames,c("Beta.Initial","sigma2.eta","sigma2.eps","u0","P0","bkout"))

cat("initialization complete...\n")
  
  
  
  # Prior specfication ------------------------------------------------------
  
  #Specify number of state differences and total number of observations
  Neta <- (map_dbl(timeDiff, length) %>% sum())
  Ntotal <-  nrow(data)
  
  
  
  W = BayesInit %>%
    map_dbl("sigma2.eta") %>%
    diag()
  
  V = BayesInit %>%
    map_dbl("sigma2.eps") %>%
    diag()
  
  m0 <- BayesInit %>%
    map_dbl("u0") 
  
  C0 <- BayesInit %>%
    map_dbl("P0") %>%
    diag()
  ###Need another way to set "TT", not ncol (y)
  prior.nu.eta <- cs + 1
  prior.Gamma.eta <- diag(1, nrow=cs)
  
  c0 <- 0.01
  d0 <- 0.01
  
  ### Beta initialize
  Beta.Initial <- BayesInit %>%
    map("Beta.Initial") %>%
    do.call("cbind", .)
  
  sigma2.beta <- 10
  
  B.star <- Beta.Initial
  B2 <- matrix(B.star, ncol = cs)
  
  
  # Transform for quicker inverse calculation in KF
  # XtX <- data %>%
  #   select(all_of(predictors)) %>%
  #   as.matrix() %>%
  #   crossprod() #hitting error here for some reason, is having "integer" variables a problem?
  #when i do this cross product manually it's working? Maybe because it's adding NACCID?
  #that's not numeric? was it adding this before???
  
  #trying this write around
  # Transform for quicker inverse calculation in KF
  XtX <- data[,predictors] %>%
      as.matrix() %>%
      crossprod()
  
  sig2beta_XtX <- sigma2.beta * XtX
  
  ex <- eigen(sig2beta_XtX, symm = TRUE)
  
  
  # Tracking obj
  Beta.Track <- array(NA, dim = c(ncol(X[[1]]), cs, numits))
  vcovWish <- array(NA, dim = c(cs, cs, numits))
  sigma2.eps.Track <- matrix(NA, nrow = numits, ncol = cs)
  
  
  It <- 0
  it <- It+1
  cat("Fitting Model...\n")
  if(!silence) pb <- txtProgressBar(min = 0, max = numits, style = 3)
  for(It in 1:numits){
    
    
    i <- 1
    y.star <- lapply(1:n, function(i){
      y[[i]]-t(X[[i]] %*% B.star) #taking out linear effects
    })
    
    
    
    ####### alpha
    bout <- lapply(1:n, function(i){
      ffbs.joint(y = y.star[[i]], V = V,W = W,m0 = m0,C0 = C0, timeDiff = timeDiff[[i]])$x
    })
    
    
    
    ####### Sigma eta
    fp <- lapply(1:n, function(i){
      apply(bout[[i]], 1, function(z)diff(z)) / sqrt( timeDiff[[i]] )
      
    }) %>%
      do.call("rbind", .) %>%
      crossprod()
    
    
    W <- cIRT::riwishart(prior.nu.eta + Neta, fp + prior.Gamma.eta)
    vcovWish[,,It] <- W
    
    ###### Sigma eps
    
    or <- lapply(1:n, function(i){
      (y.star[[i]] - bout[[i]])^2
    })  %>%
      do.call("cbind", .) %>%
      apply(1, function(x)sum((x)))
    
    
    sigma2.eps.star <- sigma2.eps.Track[It,] <- sapply(or, function(x)1/rgamma(1, ((Ntotal)/2 +c0), d0+x/2))
    V <- sigma2.eps.star * diag(cs)
    
    
    
    
    # #### Beta Post
    
    
    v.star <- lapply(1:n, function(i){
      (y[[i]]-bout[[i]])
    })
    
    B.sum <- lapply(1:n, function(i){
      (v.star[[i]] %*% X[[i]])
    }) %>%
      Reduce("+", .)
    
    B.Big <- sigma2.beta * B.sum - tcrossprod(V, Beta.Initial)
    
    Sigma.Inv <- lapply(sigma2.eps.star, function(eppps)tcrossprod(ex$vectors/(ex$values + eppps)[col(ex$vectors)], ex$vectors))
    
    BetaSim <- lapply(1:cs, function(i){
      rmvnorm(1, mean = crossprod(Sigma.Inv[[i]], B.Big[i,]), sigma = Sigma.Inv[[i]] * sigma2.eps.star[i] * sigma2.beta)
    }) %>%
      do.call("rbind", .) %>%
      t()
    
    
    B.star <- Beta.Track[,,It] <- BetaSim
    
    
    if(!silence) setTxtProgressBar(pb, It)
  }
  
  list(Beta.Track = Beta.Track, Eps.Track = sigma2.eps.Track, Eta.Track = vcovWish, Initial = BayesInit$bkout)
  
}