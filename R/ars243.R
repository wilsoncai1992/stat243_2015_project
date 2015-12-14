#### PUT ALL THE INPUT VARIABLES HERE ####

# f <- function(x) dnorm(x)
# h <- NULL

# k <- 100
# n <- 1000
# domain <- c(-Inf, Inf)

#-----------------------------------------------------------------------------------------
# TEMPORARY INPUT VALUE
h <- function(x){
  return(log(dnorm(x)))
}
h <- function(x){
  return(log(dbeta(x, 3, 2)))
}
h <- function(x){
  return(log(dexp(x, 1)))
}
k <- 4
n <- 1000
#---------------------------------------------
#This is the main engine.  It takes the log of a function (h), a required number of
#values (n), and the number of points in the initial abscissae(k).  
ars243 <- function(n, f = NULL, h = NULL, k, domain = c(-Inf, Inf)){
  source('./u_method.R')
  #-----------------------------------------------------------------------------------------
  # Checking input validity
  #-----------------------------------------------------------------------------------------
  if(is.null(f) & is.null(h)){
    stop('"f" and "h" are both missing. No default value.')
  }else{
    if(!is.null(f)){
      if ( is.expression(f) ) {
        f_exp <- f
        f <- function(x) { eval(f_exp) }
      }
      if (!is.function(f)) {
        stop( '"f" has to be either expr or function' )
      }
    }
  }
  if(k < 4){
    k <- 4
  }
  
  if(is.null(h)){
    h <- function(x){
      return(log(f(x)))
    }
  }
  
  lb <- domain[1]
  ub <- domain[2]
  #-----------------------------------------------------------------------------------------
  # Main body
  #-----------------------------------------------------------------------------------------
  #'finalValues' will be what we return.  We initiate it empty.
  finalValues <- c()
  if (all.equal(domain, c(-Inf, Inf)) == TRUE){
    abscissae.grid <- seq(-5, 5, length.out = k)
  }else{
    abscissae.grid <- seq(lb, ub, length.out = k + 2)
    abscissae.grid <- abscissae.grid[abscissae.grid > lb & abscissae.grid < ub]
  }
  
  abscissae.result <- gen.abscissae(abscissae.grid, h, lb, ub)
  #-----------------------------------------------------------------------------------------
  # Remove points where derivative is zero
  if(sum(abscissae.result[,3] == 0) > 0){
    abscissae.result <- abscissae.result[abscissae.result[,3] != 0,]
  }
  #-----------------------------------------------------------------------------------------
  ############# Renee's function of check.log.concave #############
  if(check.log.concave(abscissae.result) == FALSE){
    # expand the grid
    abscissae.grid <- seq(lb, ub, length.out = 2 * k)
    abscissae.grid <- abscissae.grid[abscissae.grid > lb & abscissae.grid < ub]
    abscissae.result <- gen.abscissae(abscissae.grid, h, lb, ub)
    
    if(check.log.concave(abscissae.result) == FALSE){
      # If still not log concave
      stop("f is not log-concave! Or k is too small.")
    }
  }
  #################################################################
  
  Tk = abscissae.result[,1]
  h_Tk = abscissae.result[,2]
  #'Coefficients' are the slope and y intercept associated with each chord in the 'l'
  #function  They are in a matrix that gets created by 'lupdater'.
  coefficients <- lupdater(Tk,h_Tk)
  # Run the code until we hit the desired length
  while(length(finalValues) < n){
    # WILSON
    Tk = abscissae.result[,1]
    h_Tk = abscissae.result[,2]
    #Take a random point from the 'sk' function, known as 'rs'.
    z = compute_z(abscissae.result)
    norm.constant = compute_norm_constant(abscissae.result,z, lb, ub)
    sampler <- rs(1, S_inv, abscissae.result, z, norm.constant, lb, ub)[1,1]
    #If that point lies outside the bounds set by my Tk values
    if(sampler < Tk[1] | sampler > Tk[length(Tk)]){
      # In this case, we need to evalue h,h'
      hValue <- h(sampler)
      h.deriv <- grad(func = h, x = sampler)
      
      # Update abscissae.result and coefficients.
      abscissae.result = update(abscissae.result,sampler,hValue,h.deriv)
      coefficients = update_coeff(coefficients,sampler,abscissae.result)
      
      ############# Renee's function of check.log.concave #############
      if(check.log.concave(abscissae.result) == FALSE){
        stop("f is not log-concave!")
      }
      #################################################################
      
    }else{ 
      #First we run a binary search to find which chord the point finds itself within.
      index <- binary(Tk, sampler)
      
      #The value of the lower bound is calculated from this chord, found with binary.
      lval <- lowerbound(sampler,coefficients,index)
      
      #The value of the upper bound is calculated by Wilson.
      uval <- du.unnormalized(sampler, abscissae.result, z, lb, ub)
      
      #We also grab a number from the uniform distribution.
      uniform <- runif(1)  
      
      #The big IF statements, this one's the 'squeezing' step
      if(uniform <= exp(lval - uval)){
        finalValues <- c(finalValues, sampler)
        
      }else{
        # Only evaluate the log of the function if we fail to squeeze.  Rejection step.
        hValue <- h(sampler)
        h.deriv <- grad(func = h, x = sampler)
        abscissae.result = update(abscissae.result,sampler,hValue,h.deriv)
        coefficients = update_coeff(coefficients,sampler,abscissae.result)
        if(uniform <= exp(hValue - uval)){
          finalValues <- c(finalValues, sampler)
          
          ############# Renee's function of check.log.concave #############
          if(check.log.concave(abscissae.result) == FALSE){
            stop("f is not log-concave!")
          }
          #################################################################
        }
      }
    }
  }
  #I was printing the length of Tk to see how many points I updated with (typically 
  #it ends up being 15-30).
  print(paste('Total number of abscissae grid point:', length(Tk)) )
  return(finalValues)
}

#-------------------------------------------------
# ars243(n = 1e3, h = h, k = 4, domain = c(0,1))
h = function(x) log(dbeta(x, 3, 2))
exp_h = function(x) dbeta(x,3,2)
curve(exp_h, from = 0, to = 1) 
lines(density(ars243(n = 1e4, h = h, k = 4, domain = c(0,1))), lty=2, col='blue')
lines(density(ars243(n = 1e4, h = h, k = 300, domain = c(0,1))), lty=2, col='red')
curve(dnorm(x), from = -2, to = 2) 
lines(density(ars243(n = 1e5, h = h, k = 4)), lty=2, col='blue')
lines(density(ars243(n = 1e4, h = h, k = 300)), lty=2, col='red')
curve(dexp(x, 1), from = 0, to = 2) 
lines(density(ars243(n = 1e4, h = h, k = 4, domain = c(0,2))), lty=2, col='blue')
lines(density(ars243(n = 1e4, h = h, k = 300, domain = c(0,2))), lty=2, col='red')

h_before_truncation = function(x){
  set.seed(0)
  # This is the data
  n=10000
  sigma=1
  theta = rnorm(n,0,sigma^2)
  total = 0
  d = 1.06*sigma*n^{-0.2}
  for(i in 1:n){
    # This is the kernal function of the data.
    total = total + dnorm((x-theta[i])/d,0,1)/(d*1000)
  }
  return(log(total))
}
h <- function(x) h_before_truncation(x)/integrate(h_before_truncation,-1,1)$value
exp_h = function(x) exp(h(x))
curve(h, from = -1, to = 1) 
# curve(exp_h, from = -1, to = 1) 
lines(density(ars243(n = 1e4, h = h, k = 3, domain = c(-1,1))), lty=2, col='blue')
lines(density(ars243(n = 1e4, h = h, k = 300, domain = c(-1,1))), lty=2, col='blue')

h = function(x) log(dnorm(x,0,1))
exp_h = function(x) dnorm(x,0,1)
curve(exp_h) 
lines(density(ars243(n = 1e4, h = h, k = 3)), lty=2, col='blue')
lines(density(ars243(n = 1e4, h = h, k = 300, domain = c(-Inf,Inf))), lty=2, col='blue')
