which.min.left <- function(x, tol = 1e-16){
  index <- which.min(x)
  if(index == 1) return(1)
  if((x[index - 1] - x[index]) < tol) return(index - 1)
  return(index)
}
