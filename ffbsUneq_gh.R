
#Forward filtering backward sampling.
ffbs.uneq2 = function(y,V,W,m0,C0, timeDiff){ 
  T = ncol(y); n = nrow(y)
  a = matrix(0, n, T); R = matrix(0, n, T) 
  m = matrix(0, n, T); C = matrix(0, n, T); B = matrix(0, n, T-1)
  H = matrix(0, n, T-1); mm = matrix(0, n, T); CC = matrix(0, n, T) 
  x = matrix(0, n, T); llike = 0.0 
  is.na.y <- is.na(y)
  t <- 0
  t <- t+1
  for (t in 1:T){ 
    if(t==1){
      a[,1] = m0; R[,1] = C0 + W
    }else{
      a[,t] = m[,t-1]; R[,t] = C[,t-1] + timeDiff[t-1, ] * W
    }
    f = a[,t]
    Q = R[,t] + V 
    A = R[,t]/Q 
    Av = ifelse(is.na(y[,t]), 0, A * (y[,t]-f))
    m[,t] = a[,t]+Av 
    C[,t] = R[,t]-Q*A**2 
    B[,t-1] = C[,t-1]/R[,t] 
    H[,t-1] = C[,t-1]-R[,t]*B[,t-1]**2 
    llike = llike + sum(dnorm(y[,t],f,sqrt(Q),log=TRUE), na.rm = TRUE) 
  }
  
  
  
  
  # C[is.na.y] <- 0
  # B[is.na.y[,-1]] <- 0
  # H[is.na.y[,-1]] <- 0
  # R[is.na.y] <- 1
  # m <- m * !is.na.y
  # a <- a * !is.na.y
  
  mm[,T] = m[,T]; CC[,T] = C[,T]
  
  C <- C * !is.na.y
  B <- B * !is.na.y[,-1]
  
  x[,T] = rnorm(n,m[,T],sqrt(C[,T])) 
  
  
  # x[,T] = m[,T]
  t <- T-1
  for (t in (T-1):1){ 
    mm[,t] = m[,t] + C[,t]/R[,t+1]*(mm[,t+1]-a[,t+1])
    CC[,t] = C[,t] - (C[,t]^2)/(R[,t+1]^2)*(R[,t+1]-CC[,t+1])
    x[,t] = rnorm(n,m[,t]+B[,t]*(x[,t+1]-a[,t+1]),sqrt(H[,t] + C[,t] * (H[,t] == 0)))
    # x[,t] = m[,t]+B[,t]*(x[,t+1]-a[,t+1])
  } 
  return(list(x=x,m=m,C=C,mm=mm,CC=CC,llike=llike)) 
}


# bk2ffbs <- list(y = y, C = C, B = B, m = m, H = H)