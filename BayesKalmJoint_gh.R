
# This program is the top level function for the multivariate state-space 
# model code. In context of the factor state-space model, this code
# is used to estimate the linear effects of covariates so that they can
# be removed from the outcomes before running the factor state-space model.

# Inputs:
# data (a data frame with all outcomes and predictors)
# outcomes (a vector of character names of outcome variables, 
# i.e. c("LOGIMEM_3","DIGIF_3","DIGIB_3"))
# predictors (a vector of character names of predictor variables/
# covariates, i.e. c("SEX_1","EDUC","NACCAGEB"))
# timevar (a character string with the column name to be used as the time,
# i.e. "VISIT_DATE")
# id (a character string with the column name to be used as the ID 
# variable, i.e. id_var="NACCID")
# initialization (type of initialization, default="Bayes")
# numits (a scalar number of iterations to perform of the Gibbs sampler)
# silence (a TRUE/FALSE value for whether or not to suppress printing
# of progress as the code runs, default=FALSE; printing not suppressed)
# seed (a scalar value to set the random seed, default=NULL)
# numitInit (a scalar number of samples to use in the Bayesian initialization
# step, default=1500)
# burnInit (a scalar number of burn in samples to discard, default=500)

#libraries and required sub-functions
library(mvtnorm)
library(tidyverse)
source("BayesKalmUneq_gh.R")
source("ffbsJoint_gh.R")

BayesKalmJoint <- function(data, outcomes, predictors, 
                                 timevar, id, initialization = "Bayes", 
                                 numits = 1000, silence = FALSE, 
                                 seed = NULL, numitInit = 1500, 
                                 burnInit = 500) {
  
  #initializing values
  cs <- length(outcomes)
  if(!is.null(seed))set.seed(seed)
  
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

    
  # Bayesian initialization
  if(initialization == "Bayes"){
    
    #running this function for each outcome one at a time
    BayesInit <- map(outcomes, function(outcomes){
      
      #creating initial model predicting outcome with all predictors
      initmodel <- lm(as.formula(paste(outcomes, "~", 
                          paste(predictors, collapse = "+"))), data = data)
      assign("initmodel",initmodel,envir=.GlobalEnv)
      
      #Storing key output from initial model
      data$resid <- resid(initmodel)
      Beta.Initial <- coef(initmodel)[-1]
      u0 <- coef(initmodel)[1]
      P0 <- data %>%
        group_by(id) %>%
        filter(time == min(time)) %>%
        .[["resid"]] %>%
        var()
      
      #Now estimating parameters for model using first 1500 samples, 
      #500 burn in, these estimates will be the initial estimates
      #in the final estimation step
      BayesKalm.Uneq(
        y.long = data[[outcomes]], X.long = as.matrix(data[predictors]),
        id = data[[id]], time = data[[timevar]],
        Burn = burnInit, Its = numitInit, 
        Beta.Initial = Beta.Initial, sigma2.beta = 10, 
        u0 = u0, P0 = P0, 
        a0 = 0.01, b0 = 0.01, c0 = 0.01, d0 = 0.01, 
        silence = FALSE
      )
      list(
        Beta.Initial <- apply(bkout$Beta[,burnInit:numitInit], 1, mean), 
        sigma2.eta <- mean(bkout$sigma2.eta[burnInit:numitInit]),
        sigma2.eps <- mean(bkout$sigma2.eps[burnInit:numitInit]),
        u0 <- u0,
        P0 <- P0,
        bkout <- list(Beta = bkout$Beta, 
                      sigma2.eta = bkout$sigma2.eta, sigma2.eps = sigma2.eps)
      )
    })
    #saving estimates from above
    assign("BayesInit",BayesInit,envir=.GlobalEnv)
  }
#Setting names of initialization parameters
BayesInit<-lapply(BayesInit,setNames,
              c("Beta.Initial","sigma2.eta","sigma2.eps","u0","P0","bkout"))

#printing that initialization step is complete
cat("initialization complete...\n")
  
  ### Prior specification ###
  
  #Specify number of state differences and total number of observations
  Neta <- (map_dbl(timeDiff, length) %>% sum())
  Ntotal <-  nrow(data)
  
  #Setting priors to be the estimates from initialization above
  #in some cases just using arbitrary priors
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
  
  prior.nu.eta <- cs + 1
  prior.Gamma.eta <- diag(1, nrow=cs)
  
  c0 <- 0.01
  d0 <- 0.01
  
  Beta.Initial <- BayesInit %>%
    map("Beta.Initial") %>%
    do.call("cbind", .)
  
  sigma2.beta <- 10
  
  B.star <- Beta.Initial
  B2 <- matrix(B.star, ncol = cs)
  
  #Calculating matrices needed in later calculations for beta posterior
  #Transform for quicker inverse calculation
  XtX <- data[,predictors] %>%
      as.matrix() %>%
      crossprod()
  
  sig2beta_XtX <- sigma2.beta * XtX 
  #why are we multiplying by prior variance for beta here??
  #it's inverse is added in the thesis
  
  ex <- eigen(sig2beta_XtX, symm = TRUE)
  
  #Creating empty tracking objcets
  Beta.Track <- array(NA, dim = c(ncol(X[[1]]), cs, numits))
  vcovWish <- array(NA, dim = c(cs, cs, numits))
  sigma2.eps.Track <- matrix(NA, nrow = numits, ncol = cs)
  
  #Initializing iteration counters
  It <- 0
  it <- It+1
  
  #printing progress updates
  cat("Fitting Model...\n")
  if(!silence) pb <- txtProgressBar(min = 0, max = numits, style = 3)
  
  #Beginning Gibbs sampler
  for(It in 1:numits){
    
    i <- 1
    #taking out linear effects from y
    y.star <- lapply(1:n, function(i){
      y[[i]]-t(X[[i]] %*% B.star)
    })
    
    #Alpha estimation using Kalman filter
    bout <- lapply(1:n, function(i){
      ffbs.joint(y = y.star[[i]], V = V,W = W,m0 = m0,C0 = C0, 
                 timeDiff = timeDiff[[i]])$x
    })
    
    #Sigma eta estimation
    fp <- lapply(1:n, function(i){
      apply(bout[[i]], 1, function(z)diff(z)) / sqrt( timeDiff[[i]] )
    }) %>% do.call("rbind", .) %>% crossprod()
    
    #Drawing random sample from posterior of Sigma eta
    W <- cIRT::riwishart(prior.nu.eta + Neta, fp + prior.Gamma.eta)
    #Storing latest Sigma eta estimate in tracking array
    vcovWish[,,It] <- W
    
    #Sigma epsilon estimation
    or <- lapply(1:n, function(i){
      (y.star[[i]] - bout[[i]])^2
    })  %>% do.call("cbind", .) %>% apply(1, function(x)sum((x)))
    
    #Drawing random sample from posterior of Sigma epsilon
    #and storing as latest Sigma epsilon estimate in tracking array
    #Doing this for each outcome from a univariate distribution
    sigma2.eps.star <- sigma2.eps.Track[It,] <- 
      sapply(or, function(x)1/rgamma(1, ((Ntotal)/2 +c0), d0+x/2))
    #Turning all of the univariate estimated into a diagonal matrix for
    #matrix estimate of Sigma epsilon
    V <- sigma2.eps.star * diag(cs)
    
    
    #Beta estimation NEED TO VERIFY THESE CALCULATIONS MATCH WITH THESIS PG 42
    
    #subtracting alpha estimates from y's
    v.star <- lapply(1:n, function(i){
      (y[[i]]-bout[[i]])
    })
    
    B.sum <- lapply(1:n, function(i){
      (v.star[[i]] %*% X[[i]])
    }) %>% Reduce("+", .)
    
    B.Big <- sigma2.beta * B.sum - tcrossprod(V, Beta.Initial)
    
    Sigma.Inv <- lapply(sigma2.eps.star, 
                        function(eppps)tcrossprod(ex$vectors/
                            (ex$values + eppps)[col(ex$vectors)], ex$vectors))
    
    #Drawing from posterior distribution of beta to get estimate
    BetaSim <- lapply(1:cs, function(i){
      rmvnorm(1, mean = crossprod(Sigma.Inv[[i]], B.Big[i,]), 
              sigma = Sigma.Inv[[i]] * sigma2.eps.star[i] * sigma2.beta)
    }) %>% do.call("rbind", .) %>% t()

    #Storing latest beta estimates in tracking matrix
    B.star <- Beta.Track[,,It] <- BetaSim
    
    #printing progress
    if(!silence) setTxtProgressBar(pb, It)
  }
  
  #storing estimates of interest
  list(Beta.Track = Beta.Track, Eps.Track = sigma2.eps.Track, 
       Eta.Track = vcovWish, Initial = BayesInit$bkout)
  
}
