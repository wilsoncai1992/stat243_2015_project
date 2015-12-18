

#---------------------------------------------
#This is the main engine.  It takes the log of a function (h), a required number of
#values (n), and the number of points in the initial abscissae(k).  
ars243 <- function(n, f = NULL, h = NULL, k, domain = c(-Inf, Inf)){
  source('./u_method.R')
  #-----------------------------------------------------------------------------------------
  # Checking input validity
  #-----------------------------------------------------------------------------------------
  M <- n/10
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
    
    while(length(finalValues) < n){
      Tk = abscissae.result[,1]
      h_Tk = abscissae.result[,2]
      coefficients <- lupdater(Tk,h_Tk)
      #Take a random point from the 'sk' function, known as 'rs'.
      z = compute_z(abscissae.result)
      norm.constant = compute_norm_constant(abscissae.result,z, lb, ub)
      # Wilson
      sampler <- as.vector(rs(M, S_inv, abscissae.result, z, norm.constant, lb, ub))
      #If that point lies outside the bounds set by my Tk values
      condition1 =  (sampler < Tk[1] | sampler > Tk[length(Tk)])
      sampler_rejected_in_1 = sampler[condition1]
      sampler_accepted_in_1 = sampler[condition1 == FALSE]
      #The value of the lower bound is calculated from this chord, found with binary.
      lval <- sapply(sampler_accepted_in_1, function(sampler) lowerbound(sampler,coefficients,binary(Tk, sampler)))
      #The value of the upper bound is calculated by Wilson.
      uval <- sapply(sampler_accepted_in_1, function(x) du.unnormalized(x, abscissae.result, z, lb, ub))
      uniform = runif(length(sampler_accepted_in_1))
      
      condition2 = (uniform <= exp(lval - uval)) 
      sampler_rejected_in_2 = sampler_accepted_in_1[condition2 == FALSE] 
      sampler_accepted_in_2 = sampler_accepted_in_1[condition2]
      uniform_rejected_in_2 = uniform[condition2 == FALSE] 
      uval_rejected_in_2 = uval[condition2 == FALSE]
      # Phase 2
      sampler_rejected_in_phase_1 = c(sampler_rejected_in_1,sampler_rejected_in_2)
      h_rejected_in_phase_1 = h(sampler_rejected_in_phase_1) 
      h.deriv <- grad(func = h, x = sampler_rejected_in_phase_1) 
      uniform_in_phase_2 = c(runif(length(sampler_rejected_in_1)),uniform_rejected_in_2)
      if(length(sampler_rejected_in_1)!=0){
        uval = c(sapply(sampler_rejected_in_1, function(x) du.unnormalized(x, abscissae.result, z, lb, ub)),
                 uval_rejected_in_2)
        
      }else{
        uval = uval_rejected_in_2
      }
      if(length(uval!=0)){
        condition2 = (uniform_in_phase_2 < exp(h_rejected_in_phase_1-uval))
        sampler_accepted_in_phase_2 = sampler_rejected_in_phase_1[condition2]
      }else{
        sampler_accepted_in_phase_2=c()
      }
      # Add accepted samples.  
      finalValues=c(finalValues,sampler_accepted_in_1,sampler_accepted_in_phase_2)
      abscissae.result = update(abscissae.result,sampler_rejected_in_phase_1,h_rejected_in_phase_1,h.deriv)
      
      if(check.log.concave(abscissae.result) == FALSE){
        stop("f is not log-concave!")
      }
    }
    #I was printing the length of Tk to see how many points I updated with (typically 
    #it ends up being 15-30).
    print(paste('Total number of abscissae grid point:', length(Tk)) )
  }
  return(finalValues)
}