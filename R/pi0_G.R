pi0_G <- function(pv, S){

  # Step 1:
  p <- length(pv)
  q <- S[1]
  gamma <- q

  tau <- rep(0, 100)
  tau0 <- gamma
  tau[1] <- tau0 + 0.5 * (0.5 - tau0)
  tau[100] <- 0.5

  d <- (tau[100] - tau[1]) / 100

  for (j in 2:99){
    tau[j] <- tau[j - 1] + d
  }

  # Step 2:

  lambda <- 1:100

  for (j in 1:100) {
    ind <- which(S <= tau[j])
    t <- S[ind]
    lambda[j] <- max(t)
  }

  beta <- 1:100

  for (j in 1:100) {
    beta[j] <- (1 / ((1 - tau[j]) * p)) + (1 / p) * (sum(pv > lambda[j]) / (1 - lambda[j]))
  }

  beta <- min(beta, 1)

  # Step 3:
  pi0 <- mean(beta)

  return(pi0)
}
