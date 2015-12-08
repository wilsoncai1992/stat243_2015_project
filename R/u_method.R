############# Renee's function of check.log.concave #############
check.log.concave <- function(abscissae.result){
  deriv.seq <- abscissae.result[,3]
  is.concave <- (deriv.seq[1] > 0) & (tail(deriv.seq, 1) < 0)
  return(is.concave)
}
#################################################################
#-----------------------------------------------------------------------------------------
# Fixed k points to create abscissae
gen.abscissae <- function(abscissae.grid, h){
  library(numDeriv)
  
  h.x <- h(abscissae.grid)
  h.deriv <- grad(func = h, x = abscissae.grid)
  abscissae.result <- cbind(abscissae.grid, h.x, h.deriv)
  return(abscissae.result)
}
#-----------------------------------------------------------------------------------------
# # Integate over exponentiated upper envelope

# Compute the norm.constant.
# compute_norm_constant = function(abscissae.result, z){
#   all.mass <- sum(
#     rep(c(-1, 1), length(abscissae.result[,3])) *
#       rep(1/abscissae.result[,3], each = 2) * 
#       exp(rep(abscissae.result[,2], each = 2) + 
#             rep(abscissae.result[,3], each = 2) * 
#             (c(-Inf, rep(z, each = 2), Inf) - rep(abscissae.result[,1], each = 2)))
#   )
#   # Normalize total density such that sum to 1
#   norm.constant <- -log(all.mass)
#   return(norm.constant)
# }

compute_norm_constant = function(abscissae.result,z, lb, ub){
  all.mass <- sum(
    rep(c(-1, 1), length(abscissae.result[,3])) *
      rep(1/abscissae.result[,3], each = 2) * 
      exp(rep(abscissae.result[,2], each = 2) + 
            rep(abscissae.result[,3], each = 2) * 
            (c(lb, rep(z, each = 2), ub) - rep(abscissae.result[,1], each = 2)))
  )
  # Normalize total density such that sum to 1
  norm.constant <- -log(all.mass)
  return(norm.constant)
}

# Compute z.
compute_z = function(abscissae.result){
  diff.h <- diff(abscissae.result[,2])
  xh.deriv <- abscissae.result[,1] * abscissae.result[,3]
  diff.xh.deriv <- diff(xh.deriv)
  diff.h.deriv <- diff(abscissae.result[,3])
  z <- (diff.h - diff.xh.deriv)/ (-diff.h.deriv)
  return(z)
}




# UN-normaized version of upper envelope (not exponentiated)
du.unnormalized <- function(x, abscissae.result, z, lb, ub){
  j <- sapply(x, function(it) min(which(it <= c(z, ub))))
  dens.u <- abscissae.result[j, 2] + (x - abscissae.result[j, 1]) * abscissae.result[j, 3]
  return(dens.u)
}

# Normaized version of upper envelope (not exponentiated).
du.normalized <- function(x, abscissae.result, z, norm.constant, lb, ub){
  j <- sapply(x, function(it) min(which(it <= c(z, ub))))
  dens.u <- abscissae.result[j, 2] + (x - abscissae.result[j, 1]) * abscissae.result[j, 3] + norm.constant
  return(dens.u)
}


#-----------------------------------------------------------------------------------------
# Inverse sampling of upper enelope density
#-----------------------------------------------------------------------------------------
# Inverse cdf of S
S_inv <- function(cdf, abscissae.result, z, norm.constant, lb, ub){
  cdf.z <- matrix(rep(c(-1, 1), length(abscissae.result[,3])) *
                    rep(1/abscissae.result[,3], each = 2) * 
                    exp(rep(abscissae.result[,2], each = 2) + 
                          rep(abscissae.result[,3], each = 2) * 
                          (c(lb, rep(z, each = 2), ub) - rep(abscissae.result[,1], each = 2)) +norm.constant),
                  byrow = TRUE, ncol = 2)
  mass.grid <- rowSums(cdf.z)
  cum.mass <- cumsum(mass.grid)
  j <- sapply(cdf, function(it) min(which(it <= cum.mass)))
  
  delta.mass <- matrix(nrow = length(j), ncol = 1)
  x.hat <- matrix(nrow = length(j), ncol = 1)
  
  z.grid <- c(z, ub)
  
  delta.mass[j==1] <- cdf[j==1]
  x.hat[j==1] <- (log(abscissae.result[1,3] * (delta.mass[j==1] - cdf.z[1,1])) - abscissae.result[1,2] - norm.constant)/abscissae.result[1,3] + abscissae.result[1,1]
  x.hat[j==1][x.hat[j==1] < lb] <- lb
  
  j.nonzero <- j[j!=1]
  delta.mass[j!=1] <- cdf[j!=1] - cum.mass[j.nonzero - 1]#[j.nonzero - 1]
  partial.formula <- exp(abscissae.result[j.nonzero,2] + abscissae.result[j.nonzero,3] * (z.grid[j.nonzero - 1] - abscissae.result[j.nonzero,1]) + norm.constant) + delta.mass[j!=1] * abscissae.result[j.nonzero,3]
  x.hat[j!=1] <- (log(partial.formula) - abscissae.result[j.nonzero,2] - norm.constant) / abscissae.result[j.nonzero,3] + abscissae.result[j.nonzero,1]
  
  j.last <- j[j==length(cum.mass)]
  # x.hat[j==length(cum.mass)] <- (log(abscissae.result[j.last,3] * (delta.mass[j==length(cum.mass)] - tail(mass.grid, 1))) - abscissae.result[j.last,2] - norm.constant)/abscissae.result[j.last,3] + abscissae.result[j.last,1]
  x.hat[j==length(cum.mass)] <- (log(abscissae.result[j.last,3] * (cdf[j==length(cum.mass)] - 1.0 + tail(as.vector(cdf.z),1))*as.integer((cdf[j==length(cum.mass)] - 1.0) <= 0)) - abscissae.result[j.last,2] - norm.constant)/abscissae.result[j.last,3] + abscissae.result[j.last,1]
  
  return(x.hat)
}


# Inverse sampling function
rs <- function(n.sim, S_inv, abscissae.result, z, norm.constant, lb, ub){
  u.temp <- runif(n = n.sim)
  x.temp <- S_inv(u.temp, abscissae.result, z, norm.constant, lb, ub)
  return(x.temp)
}


#-------------------------------------------------------------------------------------------

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
  k=length(Tk)
  coefficients1 <- matrix(c((h_Tk[1:(k-1)]-h_Tk[2:k])/(Tk[1:(k-1)]-Tk[2:k])), k-1, 2, byrow = FALSE)
  coefficients1[,2] = h_Tk[1:(k-1)] - coefficients1[,1]*Tk[1:(k-1)]
  return(coefficients-coefficients1)
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
#     coefficients[i_added-1,1] <- (h_Tk[i_added-1]-h_Tk[i_added]) / (Tk[i_added-1]-Tk[i_added])   
#     coefficients[i_added-1,2] <- h_Tk[i_added-1] - coefficients[i_added-1,1] * Tk[i_added-1]
#     new_x = (h_star-h_Tk[i_added+1]) / (x_star-Tk[i_added+1])   
#     new_row = c(new_x,h_star-new_x*x_star)
    # updated_coefficients = rbind(coefficients[1:(i_added-1),],new_row,
                                 # coefficients[i_added:length(coefficients[,1]),])
    new_x = (h_star-h_Tk[i_added+1]) / (x_star-Tk[i_added+1])   
    new_row = c(new_x,h_star-new_x*x_star)
    if(i_added<=length(coefficients[,1])){
      updated_coefficients = rbind(coefficients[1:(i_added-1),],new_row,
                                   coefficients[i_added:length(coefficients[,1]),])
    }else{
      updated_coefficients = rbind(coefficients[1:(i_added-1),],new_row)
    }
    updated_coefficients[i_added-1,1] <- (h_Tk[i_added-1]-h_Tk[i_added]) / (Tk[i_added-1]-Tk[i_added])   
    updated_coefficients[i_added-1,2] <- h_Tk[i_added-1] - updated_coefficients[i_added-1,1] * Tk[i_added-1]
    
#     updated_coefficients = rbind(coefficients[1:(i_added-2),],new_row,
#                                  coefficients[(i_added-1):length(coefficients[,1]),])
    
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
