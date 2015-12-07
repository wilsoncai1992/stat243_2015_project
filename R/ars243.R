#### PUT ALL THE INPUT VARIABLES HERE ####

# f <- dnorm()
# h <- NULL

# k <- 100
# n <- 1000
# domain <- c(-Inf, Inf)
#-----------------------------------------------------------------------------------------
# TEMPORARY INPUT VALUE
h <- function(x){
  return(log(dnorm(x)))
}
k <- 4
n <- 1000
#-----------------------------------------------------------------------------------------
#### START OF FUNCTION ####
source('./u_method.R')
if(is.null(h)){
  h <- function(x){
    return(log(f(x)))
  }
}

#---------------------------------------------
#This is the main engine.  It takes the log of a function (h), a required number of
#values (n), and the number of points in the initial abscissae(k).  
ars243 <- function(h, n, k){
  #'finalValues' will be what we return.  We initiate it empty.
  finalValues <- c()
  abscissae.grid <- seq(-5, 5, length.out = k)
  abscissae.result <- gen.abscissae(abscissae.grid, h)
  Tk = abscissae.result[,1]
  h_Tk = abscissae.result[,2]
  #'Coefficients' are the slope and y intercept associated with each chord in the 'l'
  #function  They are in a matrix that gets created by 'lupdater'.
  coefficients <- lupdater(Tk,h_Tk)
  # Run the code until we hit the desired length
  while(length(finalValues) < n){
    #Take a random point from the 'sk' function, known as 'rs'.
    z = compute_z(abscissae.result)
    norm.constant = compute_norm_constant(abscissae.result,z)
    sampler <- (rs(1, S_inv, abscissae.result, z, norm.constant))[1,1]
    #If that point lies outside the bounds set by my Tk values
    if(sampler < Tk[1] | sampler > Tk[length(Tk)]){
      # In this case, we need to evalue h,h'
      hValue <- h(sampler)
      h.deriv <- grad(func = h, x = sampler)
      
      # Update abscissae.result and coefficients.
      abscissae.result = update(abscissae.result,sampler,hValue,h.deriv)
      coefficients = update_coeff(coefficients,sampler,abscissae.result)
    }else{ 
      #First we run a binary search to find which chord the point finds itself within.
      index <- binary(Tk, sampler)
      
      #The value of the lower bound is calculated from this chord, found with binary.
      lval <- lowerbound(sampler,coefficients,index)
      
      #The value of the upper bound is calculated by Wilson.
      uval <- du.unnormalized(sampler, abscissae.result, z)
      
      #We also grab a number from the uniform distribution.
      uniform <- runif(1)  
      
      #The big IF statements, this one's the 'squeezing' step
      if(uniform <= exp(lval - uval)){
        finalValues <- c(finalValues, sampler)
        
      }else{
        # Only evaluate the log of the function if we fail to squeeze.  Rejection step.
        hValue <- h(sampler)
        h.deriv <- grad(func = h, x = sampler)
        if(uniform <= exp(hValue - uval)){
          finalValues <- c(finalValues, sampler)
          abscissae.result = update(abscissae.result,sampler,hValue,h.deriv)
          coefficients = update_coeff(coefficients,sampler,abscissae.result)
        }
      }
    }
  }
  #I was printing the length of Tk to see how many points I updated with (typically 
  #it ends up being 15-30).
  print(length(coefficients[,1]))
  return(finalValues)
}

#-------------------------------------------------
ars243(h, n = 1e3, k)
