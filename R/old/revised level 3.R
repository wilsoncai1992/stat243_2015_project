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
k <- 100
n <- 1000
#-----------------------------------------------------------------------------------------
#### START OF FUNCTION ####
if(is.null(h)){
  h <- function(x){
    return(log(f(x)))
  }
}
#-----------------------------------------------------------------------------------------
#
# Fixed k points to create abscissae
abscissae.grid <- seq(-5, 5, length.out = k)
gen.abscissae <- function(abscissae.grid, h){
  library(numDeriv)
  h.x <- h(abscissae.grid)
  h.deriv <- grad(func = h, x = abscissae.grid)
  abscissae.result <- cbind(abscissae.grid, h.x, h.deriv)
  return(abscissae.result)
}
diff.h <- diff(abscissae.result[,2])
xh.deriv <- abscissae.result[,1] * abscissae.result[,3]
diff.xh.deriv <- diff(xh.deriv)
diff.h.deriv <- diff(abscissae.result[,3])

z <- (diff.h - diff.xh.deriv)/ (-diff.h.deriv)
#-----------------------------------------------------------------------------------------
# Integate over exponentiated upper envelope
all.mass <- sum(
  rep(c(-1, 1), length(abscissae.result[,3])) *
    rep(1/abscissae.result[,3], each = 2) * 
    exp(rep(abscissae.result[,2], each = 2) + 
          rep(abscissae.result[,3], each = 2) * 
          (c(-Inf, rep(z, each = 2), Inf) - rep(abscissae.result[,1], each = 2)))
)
# Normalize total density such that sum to 1
norm.constant <- -log(all.mass)

# UN-normaized version of upper envelope (not exponentiated)
du.unnormalized <- function(x, abscissae.result, z){
  j <- sapply(x, function(it) min(which(it <= c(z, Inf))))
  dens.u <- abscissae.result[j, 2] + (x - abscissae.result[j, 1]) * abscissae.result[j, 3]
  return(dens.u)
}

# Normaized version of upper envelope (not exponentiated).
du.normalized <- function(x, abscissae.result, z, norm.constant){
  j <- sapply(x, function(it) min(which(it <= c(z, Inf))))
  dens.u <- abscissae.result[j, 2] + (x - abscissae.result[j, 1]) * abscissae.result[j, 3] + norm.constant
  return(dens.u)
}


#-----------------------------------------------------------------------------------------
# Inverse sampling of upper enelope density
#-----------------------------------------------------------------------------------------
# Inverse cdf of S
S_inv <- function(cdf, abscissae.result, z, norm.constant){
  mass.grid <- rowSums(matrix(rep(c(-1, 1), length(abscissae.result[,3])) *
                                rep(1/abscissae.result[,3], each = 2) * 
                                exp(rep(abscissae.result[,2], each = 2) + 
                                      rep(abscissae.result[,3], each = 2) * 
                                      (c(-Inf, rep(z, each = 2), Inf) - rep(abscissae.result[,1], each = 2)) +norm.constant),
                              byrow = TRUE, ncol = 2))
  cum.mass <- cumsum(mass.grid)
  j <- sapply(cdf, function(it) min(which(it <= cum.mass)))
  
  delta.mass <- matrix(nrow = length(j), ncol = 1)
  x.hat <- matrix(nrow = length(j), ncol = 1)
  
  z.grid <- c(z, Inf)
  
  delta.mass[j==1] <- cdf[j==1]
  x.hat[j==1] <- (log(abscissae.result[1,3] * delta.mass[j==1]) - abscissae.result[1,2] - norm.constant)/abscissae.result[1,3] + abscissae.result[1,1]
  
  j.nonzero <- j[j!=1]
  delta.mass[j!=1] <- cdf[j!=1] - cum.mass[j.nonzero - 1]#[j.nonzero - 1]
  partial.formula <- exp(abscissae.result[j.nonzero,2] + abscissae.result[j.nonzero,3] * (z.grid[j.nonzero - 1] - abscissae.result[j.nonzero,1]) + norm.constant) + delta.mass[j!=1] * abscissae.result[j.nonzero,3]
  x.hat[j!=1] <- (log(partial.formula) - abscissae.result[j.nonzero,2] - norm.constant) / abscissae.result[j.nonzero,3] + abscissae.result[j.nonzero,1]
  
  j.last <- j[j==length(cum.mass)]
  # x.hat[j==length(cum.mass)] <- (log(abscissae.result[j.last,3] * (delta.mass[j==length(cum.mass)] - tail(mass.grid, 1))) - abscissae.result[j.last,2] - norm.constant)/abscissae.result[j.last,3] + abscissae.result[j.last,1]
  x.hat[j==length(cum.mass)] <- (log(abscissae.result[j.last,3] * (cdf[j==length(cum.mass)] - 1.0)*as.integer((cdf[j==length(cum.mass)] - 1.0) <= 0)) - abscissae.result[j.last,2] - norm.constant)/abscissae.result[j.last,3] + abscissae.result[j.last,1]
  
  return(x.hat)
}


# Inverse sampling function
rs <- function(n.sim, S_inv, abscissae.result, z, norm.constant){
  u.temp <- runif(n = n.sim)
  x.temp <- S_inv(u.temp, abscissae.result, z, norm.constant)
  return(x.temp)
}

#-----------------------------------------------------------------------------------------
# rs(n.sim, S_inv, abscissae.result, z, norm.constant)
x <- seq(-6, 6, by = 0.1)
du.normalized(x, abscissae.result, z, norm.constant)

haha2 <- function(x){ exp(du.unnormalized(x, abscissae.result, z))}
curve(haha2, from = -3, to = 3)
haha <- function(x){ exp(du.normalized(x, abscissae.result, z, norm.constant))}
curve(haha, from = -3, to = 3)

lines(density(rs(n, S_inv, abscissae.result, z, norm.constant)), lty = 2 , col = 'red')
lines(density(rs(1e5, S_inv, abscissae.result, z, norm.constant)), lty = 2 , col = 'red')









# Set up arguments:
  n <- 1000

 
#---------------------------------------------

#Example function for f
h <- function(x){
  return(log(dnorm(x, 0, 1)))}
n=10
k=4

#---------------------------------------------
#This is the main engine.  It takes the log of a function (h), a required number of
#values (n), and the number of points in the initial abscissae(k).  
level3 <- function(h, n, k){
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
  print(length(Tk))
  return(finalValues)
}

#----------------------------------------------
#This function runs a binary search to find the chord we should use
#to represent this point within the lower function. 

binary <- function(Tk, sampledPoint){
  bot <- 1
  top <- (length(Tk))
  for(k in 1:(log2(length(Tk)))){
    q <- Tk[(bot+(top-bot)/2)]
    if(sampledPoint > q){
      bot <- (bot+(top-bot)/2)
    }else{
      top <- (bot+(top-bot)/2)}
  }
  index <- floor(((top+bot)/2))
  return(index)
}

#-------------------------------------------------

#This function runs the 'l' function.
lowerbound <- function(x, coefficients, index){
  lvalue <- coefficients[index,1]*x + coefficients[index,2]
  return(lvalue)
}
#-------------------------------------------------
# This function updates the chords of the 'l' function with new Tk points.
# The inputs are (sorted) Tk and the values of h at points of Tk respectively, 
# represented as two vectors. 
lupdater <- function(Tk,h_Tk){
  coefficients <- matrix(0, nrow = (length(Tk)-1), ncol = 2)
  for(j in 1:(length(Tk)-1)){
    coefficients[j,1] <- (h_Tk[j]-h_Tk[j+1])/(Tk[j]-Tk[j+1])   
    coefficients[j,2] <- h_Tk[j] - coefficients[j,1]*Tk[j]
  }
  return(coefficients)
}

# This function performs the updating step.
# The inputs are original_abs, which is a matrix of dimension k times 3, storing x,h,h'..
# and the original coefficients matrix.
# and the point to be added, x_star, and h_star=h(x_star), h_p_star = h'(x_star).
# The outputs are the coordinates and derivative of abscissae with a new point added, as a matrix with three columns, each recording x_coord, y_coord and derivatives.
update <- function(original_abs,x_star,h_star,h_deri_star){
  # Represent the information about x_star as a vector.
  added_abs = c(x_star,h_star,h_deri_star)
  # Add the new point to the original abscissae.
  updated_abs = rbind(added_abs,original_abs)
  # Sort by the x coordinates in ascending order.
  updated_abs = updated_abs[order(updated_abs[,1]),]
  return(updated_abs)
}

# This function updates the coefficient matrix.
update_coeff <- function(coefficients,x_star,updated_abscissae.result){
  # The index of the added point:
  Tk = updated_abscissae.result[,1] 
  h_Tk = updated_abscissae.result[,2] 
  i_added = which(Tk == x_star)
  h_star = h_Tk[i_added]
  
  # Update the coefficient matrix.
  if(i_added<length(Tk) && i_added>1){
    coefficients[i_added-1,1] <- (h_Tk[i_added-1]-h_Tk[i_added]) / (Tk[i_added-1]-Tk[i_added])   
    coefficients[i_added-1,2] <- h_Tk[i_added-1] - coefficients[i_added-1,1] * Tk[i_added-1]
    new_x = (h_star-h_Tk[i_added+1]) / (x_star-Tk[i_added+1])   
    new_row = c(new_x,h_star-new_x*x_star)
    updated_coefficients = rbind(coefficients[1:(i_added-1),],new_row,
                                 coefficients[i_added:length(coefficients[,1]),])
    
  }
  
  if(i_added == 1){
    new_x = (h_star-h_Tk[i_added+1]) / (x_star-Tk[i_added+1])   
    new_column = c(new_x,h_star-new_x*x_star)
    updated_coefficients = rbind(new_column, coefficients)
  }
  if(i_added == length(Tk)){
    new_x = (h_star-h_Tk[length(h_Tk)-1]) / (x_star-Tk[length(Tk)-1])   
    new_column = c(new_x,h_star-new_x*x_star)
    updated_coefficients = rbind(coefficients,new_column)
  }
  
  return(updated_coefficients)
}



# Test the update_coeffi function.
k=4
finalValues <- c()
abscissae.grid <- seq(-5, 5, length.out = k)
abscissae.result <- gen.abscissae(abscissae.grid, h)
Tk = abscissae.result[,1]
h_Tk = abscissae.result[,2]
#'Coefficients' are the slope and y intercept associated with each chord in the 'l'
#function  They are in a matrix that gets created by 'lupdater'.
coefficients <- lupdater(Tk,h_Tk)
abscissae.result = update(abscissae.result,sampler,hValue,h.deriv)
coefficients = update_coeff(coefficients,x_star=sampler,updated_abscissae.result=abscissae.result)
# This function packs Wilson's code to compute z.
compute_z = function(abscissae.result){
  diff.h <- diff(abscissae.result[,2])
  xh.deriv <- abscissae.result[,1] * abscissae.result[,3]
  diff.xh.deriv <- diff(xh.deriv)
  diff.h.deriv <- diff(abscissae.result[,3])
  z <- (diff.h - diff.xh.deriv)/ (-diff.h.deriv)
  return(z)
}

# This function uses Wilson's code to compute the norm.constant.
compute_norm_constant = function(abscissae.result,z){
  all.mass <- sum(
    rep(c(-1, 1), length(abscissae.result[,3])) *
      rep(1/abscissae.result[,3], each = 2) * 
      exp(rep(abscissae.result[,2], each = 2) + 
            rep(abscissae.result[,3], each = 2) * 
            (c(-Inf, rep(z, each = 2), Inf) - rep(abscissae.result[,1], each = 2)))
  )
  # Normalize total density such that sum to 1
  norm.constant <- -log(all.mass)
  return(norm.constant)
}

level3(h,10,k)