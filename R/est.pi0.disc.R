est.pi0.disc <- function(p.count, p.grid, n = 20, lambda.max = 0.5){

  m <- sum(p.count)

  cum.count <- cumsum(p.count)

  target <- 1/n*(1:(n-1))

  candidate <- p.grid[p.grid > 0 & p.grid < 1]

  # For each target, find the closest candidate
  lambda.vec <- unique(sapply(target, function(x) candidate[which.min.left(abs(candidate - x))]))
  if(length(lambda.vec) == 1) return(list(pi0 = (m - cum.count[which(lambda.vec == p.grid)]) / (1 - lambda.vec)/m, lambda = lambda.vec))
  if(lambda.max < 1) lambda.vec <- lambda.vec[lambda.vec <= lambda.max]
  gap <- c(lambda.vec[1], diff(lambda.vec))

  index.lambda <- sapply(lambda.vec, function(x) which(p.grid == x))
  cum.bin.count <- cum.count[index.lambda]

  # All but the last bin, length B-1
  bin.counts <- c(cum.bin.count[1], diff(cum.bin.count))

  # Number of bins
  B <- length(lambda.vec)+1
  R <- cumsum(bin.counts)
  tail.m0 <- (m - R) / (1 - lambda.vec)
  temp <- bin.counts/gap - tail.m0

  if(sum(temp <= 0) > 0){
    index <- min((1:(B - 1))[temp <= 0])
  }else{
    index <- B-1
  }
  return(list(pi0 = tail.m0[index] / m, lambda = lambda.vec[index]))

}
