q.value_pi0 <- function(pv, pi0){
  p <- length(pv)
  qvalues <- 1:p
  pv_s <- sort(pv)
  ind <- which(pv == pv_s[p])
  qvalues[p] <- pi0 * pv_s[p]
  qvalues_S <- rep(0, p)
  qvalues_S[ind] <- qvalues[p]
  seq <- rev(1:(p - 1))

  for (i in (seq)) {
    ind <- which(pv == pv_s[i])
    a <- (pi0 * p * pv_s[i]) / i
    qvalues[i] <- min(qvalues[i + 1], a)
    qvalues_S[ind] <- qvalues[i]
  }
  return(qvalues_S)
}
