get.tbl.full <- function(x, NN, ss){
  tbl <- table(x)
  if(length(tbl) == NN) return(tbl)
  res <- rep(0, NN)
  names(res) <- ss
  res[names(tbl)] <- tbl
  return(res)
}
