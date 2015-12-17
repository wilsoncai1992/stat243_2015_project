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
  # If uniform
  if(all(abscissae.result[,3] == 0)){
    finalValues <- runif(n, min = lb, max = ub)
  #-----------------------------------------------------------------------------------------
  }else{
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
        # coefficients = update_coeff(coefficients,sampler,abscissae.result)
        Tk = abscissae.result[,1]
        h_Tk = abscissae.result[,2]
        coefficients <- lupdater(Tk,h_Tk)
        ############# Renee's function of check.log.concave #############
        if(check.log.concave(abscissae.result) == FALSE){
          stop("f is not log-concave!")
        }
        #################################################################
        #First we run a binary search to find which chord the point finds itself within.
        index <- binary(Tk, sampler)
        
        #The value of the lower bound is calculated from this chord, found with binary.
        lval <- lowerbound(sampler,coefficients,index)
        
        #The value of the upper bound is calculated by Wilson.
        uval <- du.unnormalized(sampler, abscissae.result, z, lb, ub)
        
        #We also grab a number from the uniform distribution.
        uniform <- runif(1)  
        
        hValue <- h(sampler)
        h.deriv <- grad(func = h, x = sampler)
        if(uniform <= exp(hValue - uval)){
          finalValues <- c(finalValues, sampler)
        }
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
          # coefficients = update_coeff(coefficients,sampler,abscissae.result)
          Tk = abscissae.result[,1]
          h_Tk = abscissae.result[,2]
          coefficients <- lupdater(Tk,h_Tk)
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
  }
  return(finalValues)
}