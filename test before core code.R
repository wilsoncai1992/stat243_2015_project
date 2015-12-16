### derivative at left end point should be greater or equal to that of the right end point.
log_concavity_check <- function(abscissae.result){
  n=length(abscissae.result[,3])
  return(abscissae.result[1,3]-abscissae.result[n,3]>=-1e-8)
}


