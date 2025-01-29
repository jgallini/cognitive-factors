
#This program creates a function that estimates factor scores for each 
#participant at each time point (alphas), the factor loading matrix G, 
#and the variance components sigma eta and sigma epsilon

# Inputs to this function: 
# 1. data (a data frame), 
# 2. outcome_vars (a vector of character names of outcome variables, 
# i.e. c("LOGIMEM_3","DIGIF_3","DIGIB_3")
# 3. id_var (a character string with the column name to be used as the ID 
#variable, i.e. id_var="NACCID")
# 4. time_var (a character string with the column name to be used as the 
# time variable, i.e. time_var="VISIT_DATE")
# 5. G-form (a matrix with 0's and 1's specifying which of the outcome 
# variables can load onto which of the factors)
# 6. sigma2.eps.star (a vector of prior point estimates for 
# sigma squared epsilon, should be 1 x number of outcome variables. Can be 
# derived from estimates from the multivariate model.)
# 7. sigeps.allests (a matrix (dimension: iterations in multivariate model
# default 1000 x outcome_vars) of all estimates for sigma epsilon squared
# from the multivariate run of the model, used for prior parameter 
# calculations for sigma epsilon squared distributions)
# 8. its (a numeric value for the number of iterations to run, default 10000)
# 9. burn (numeric scalar value of samples to burn, default floor(its/2))
# 10. wait (numeric scalar value of number of samples to wait before beginning 
# variance calculations, default floor(its/4))
# 11. mu_G (a matrix of prior means for G, must be of dimension 
# outcome_vars x number of factors, defaults to matrix of 0's)
# 12. var_G (a numeric scalar value for the prior variance of G elements)
# 13. m0 (a numeric scalar value for the prior means of alpha, default 0)
# 14. C0 (a numeric scalar value for the prior variances of alpha, default 10)
# 15. seed (a numeric scalar value to set the random seed, default 5631)

# Function outputs: A list with estimates from all iterations for G, 
# alpha, sigma epsilon squared, and sigma eta squared.

#calling ffbs.fac function, needed for Kalman filter calculations
library(here)
here::i_am("Code/SSM_fac_gh.R")
source(here("Code","ffbsFac_gh.R"))


SSM_fac <- function(
    data, outcome_vars, id_var, time_var, G_form, sigma2.eps.star, 
    sigeps.allests, its = 10000, burn = floor(its/2), wait = floor(its/4),
    mu_G = matrix(0, length(outcome_vars), ncol(G_form)), 
    var_G = 10, #mean and variance priors for G
    m0 = 0, C0 = 10, #priors for alpha
    seed = 5631
){
  #libraries
  library(tidyverse)
  
  #setting seed if needed
  if(!is.null(seed))set.seed(seed)
  
  #creating dimension variables to be used throughout
  K <- length(outcome_vars)
  Q <- ncol(G_form)
  
  #Initialize y,time, and alpha lists
  y <- list()
  time <- list()
  alpha.star <- list()
  
  #creating id-based parameters to use later
  uid <- unique(data[[id_var]]) #all unique id's
  nid <- N <- length(uid) #number of unique participants
  Ntotal <- nrow(data) #number of total observations
  Neta <- nrow(data) - nid #diff between total observations and number of 
  #participants (N(J-1)), will be component in posterior
  #degrees of freedom calculation for sigma eta distribution
  
  #creating matrix of outcomes,matrix of time points, and initial matrix
  #of alpha values for each individual participant, 
  #storing all matrices for all participants in respective lists
  for(i in 1:nid){
    tfid <- data[id_var] == uid[i]
    y[[i]] <- t(as.matrix(data[tfid, outcome_vars]))
    time[[i]] <- data[[time_var]][tfid]
    
    alpha.star[[i]] <- matrix(m0, nrow = nrow(y[[i]]), ncol = ncol(y[[i]]))
  }
  
  #calculating time difference variable
  time_diff <- map(time, diff)

  # Variance Priors
  #Sigma epsilon variance matrix referred to as V throughout
  V <- V_0 <- diag(sigma2.eps.star, K, K) 
  #Sigma eta variance matrix referred to as W throughout
  W <- W_0 <- diag(1,Q,Q)
  
  
  #Calculating 2 prior parameters c0 and d0 for inverse gamma distribution 
  #for each sigma squared epsilon based on point estimates from multivariate
  #model run
  
  #taking variance of each element of sigma squared epsilon across all
  #iterations from multivariate model
  var.sigeps.priors<-apply(sigeps.allests,2,var)
  #initiating vectors for prior parameters c0/see0 and d0
  see0<-vector()
  d0<-vector()
  #transforming the mean and variance of an inverse gamma back into the two
  #parameters that define an inverse gamma distribution (see0 and d0)
  for (k in 1:10){
    see0[k]<-(2*((sigma2.eps.star[k]**2)/var.sigeps.priors[k]))+4
    d0[k]<-(2*((sigma2.eps.star[k]**3)/var.sigeps.priors[k]))+
      2*sigma2.eps.star[k]
  }

  #Variance holders, empty for now 
  vcovWish <- matrix(NA, its, K)
  wcovWish <- array(NA, dim = c(Q, Q, its))
  
  #Prior degrees of freedom parameter for sigma eta distribution
  prior.nu.eta <- Q + 1
  
  #starting G estimates with the prior, G start will update throughout sampler
  G.star <- mu_G 
  
  #initiating G array for tracking all iterations of sampler
  G.track <- array(dim = c(dim(G_form), its))
  
  #initiating list for tracking cross product of alpha values
  ata.inv <- list()
  
  #creating matrix of all cognitive testing outcome data
  ybig <- t(as.matrix(data[outcome_vars]))
  
  #initiating list to keep track of alphas for each iteration
  alpha.track <- list()
  
  #beginning Gibbs sampler
  for(i in 1:its){
    
    #Estimating G transformed Alphas- running Kalman filter/smoother
    alpha.star <- alpha.track[[i]] <- lapply(1:nid, function(yi){
      ffbs.fac(
        y = y[[yi]], G = G.star, 
        V = V, W = W, 
        m0 = 0, C0 = diag(10, Q, Q), 
        timeDiff = time_diff[[yi]]
      )$x
    })
    #keeping track of where the simulation is
    print(i)
    
    
    #Estimating G
    abig <- do.call("cbind", alpha.star) #all current alpha values
    aybig <- tcrossprod(ybig, abig) #sum of y's by alpha's
    
    #solving for mean of the posterior for each k
    mu.g<-array(NA,c(K,Q))
    for (k in 1:K){
      ata.inv[[k]] <-solve(var_G*tcrossprod(abig) + 
                             diag(sigma2.eps.star[k],Q,Q))
      mu.g[k,] <- (sigma2.eps.star[k]%*%mu_G[k,] + 
                     var_G%*%aybig[k,]) %*% ata.inv[[k]]
    }
    
    #drawing a random sample from the posterior distribution of G to be
    #current G estimate
      G.track[,,i] <- G.star <- abs((lapply(1:K, function(k){
        mvtnorm::rmvnorm(1, mean = mu.g[k,], 
                         sigma = ata.inv[[k]] * sigma2.eps.star[k] * var_G)
       }) %>% do.call("rbind", .)) * G_form) 
  

    #Estimating variance components
      
    #Sigma eta  
    if(i > wait){
      
      #calculating likelihood portion of posterior parameter
      fp <- lapply(1:N, function(i){
        apply(alpha.star[[i]], 1, function(z)diff(z))/sqrt(time_diff[[i]] )
      }) %>% do.call("rbind", .) %>% crossprod()
      
      #combining with prior to form full posterior parameter
      fp2<-cov2cor(fp+W_0)
      
      #drawing a random sample from the posterior distribution of Sigma Eta (W)
      #to be the current Sigma Eta estimate
      W <- cov2cor(cIRT::riwishart(Neta+prior.nu.eta, fp2))
      
      wcovWish[,,i] <- W    
    }
    
    #Sigma epsilon
    
    #calculating likelihood portion of posterior
    or <- lapply(1:N, function(i){
      (y[[i]] - G.star %*% alpha.star[[i]])^2
    })  %>% do.call("cbind", .) %>% apply(1, function(x)sum((x))) 
    
    #for each of k outcomes, drawing a random sample from the (univariate) 
    #posterior distribution for sigma epsilon k to be current sigma epsilon
    #k estimate
    for (k in 1:10){
      sigma2.eps.star[k] <- 1/rgamma(1, (((Ntotal) + see0[k])/2), 
                                     (or[k]+d0[k])/2)
    }
    
    #storing all estimates
    vcovWish[i,] <- sigma2.eps.star
    #turning into a diagonal matrix
    V <- vcovWish[i,]* diag(K)
    
    #keeping track of where the simulation is
    if((i %% 10) == 0)print(i)
    
  } 
  
  #output from this function: G, alpha, Sigma eta, and Sigma Epsilon estimates
  list(G = G.track, alpha = alpha.track, Eps = vcovWish,
       Eta = wcovWish)
  
}
