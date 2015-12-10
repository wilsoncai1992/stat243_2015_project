### log concavity check

du.unnormalized <- function(x, abscissae.result, z, lb, ub){
  j <- sapply(x, function(it) min(which(it <= c(z, ub))))
  dens.u <- abscissae.result[j, 2] + (x - abscissae.result[j, 1]) * abscissae.result[j, 3]
  return(dens.u)
}

lowerbound <- function(x, coefficients, index){
  lvalue <- coefficients[index,1]*x + coefficients[index,2]
  return(lvalue)
}

### for upper envelop and lower envelop, the derivatives should be decreaseing
log_concavity_check <- function(f,abscissae.grid){
  abscissae.grid = uniq(sort(abscissae.grid))
  n = length(abscissae.grid)
  eps = min(abscissae.grid[2:n] - abscissae.grid[1:(n-1)])/100
  x_right = abscissae.grid + eps
  x_left = abscissae.grid + eps
  library(numDeriv)
  # derivative of upper envelop at each x
  du = grad(du.unnormalized, abscissae.grid, method = "Richardson")
  # derivative of lower envelop at left of each x
  l = grad(lowerbound, x_left, method = "Richardson")
  # derivative of lower envelop at right of each x
  r = grad(lowerbound, x_right, method = "Richardson")
  ## the derivatives should be decreasing
  result = (l[2:n] > du[2:n]) * (du[1:(n-1)] > r[1:(n-1)])
  return(prod(result))

}


