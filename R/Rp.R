
Rp <- function(pv, s){
  l_pv <- length(pv)
  pi0 <- 1:500

  for (i in 1:500){
    pv_R <- 1:l_pv
    for(j in 1:l_pv){

      u <- runif(1)
      eso <- rep(pv[j], length(s))
      esa <- abs(s - eso)
      min_esa <- min(esa)
      ind <- which(esa == min_esa)
      if (ind > 1) {
        pv_R[j]<-pv[j]-u*(pv[j]-s[ind-1])
        if(pv_R[j]<0) {pv_R[j]<-0
        }
        }
          if (ind == 1) {
            pv_R[j] <- pv[j] - u * (pv[j] - 0)
            if(pv_R[j]<0) {pv_R[j]<-0
            }
          }

    }
    lambda <- 0.5
    pi0[i] <- sum(pv_R > lambda) / ((1 - lambda) * length(pv_R))
  }
  return(mean(pi0))
}
