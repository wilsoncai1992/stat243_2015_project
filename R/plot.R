
set.seed(0)
k=300
h = function(x) log(dbeta(x, 3, 2))
abscissae.grid <- seq(lb, ub, length.out = k)
abscissae.grid <- abscissae.grid[abscissae.grid > lb & abscissae.grid < ub]
abscissae.result <- gen.abscissae(abscissae.grid, h)
Tk = abscissae.result[,1]
h_Tk = abscissae.result[,2]
coefficients <- lupdater(Tk,h_Tk)
finalValues = c()
lb <- 0
ub <- 1

#The value of the lower bound is calculated from this chord, found with binary.
lval <- function(sampler) lowerbound(sampler,coefficients,index)

#The value of the upper bound is calculated by Wilson.
uval <- function(sampler) du.unnormalized(sampler, abscissae.result, z, lb, ub)

first_step <- function(x) exp(lval(x) - uval(x))
second_step <- function(x) exp(h(x) - uval(x))

library(ggplot2)
new_point = data.frame(x=c(0.3,0.5),y=c(0.5,1))
ggplot(data.frame(x=c(0,1)), aes(x)) +
  stat_function(fun=first_step, geom="line", aes(colour="phase 1")) +
  stat_function(fun=second_step, geom="line", aes(colour="phase 2"))  
 




z = compute_z(abscissae.result)
norm.constant = compute_norm_constant(abscissae.result,z, lb, ub)
sampler <- rs(1, S_inv, abscissae.result, z, norm.constant, lb, ub)[1,1]
index <- binary(Tk, sampler)
#The value of the lower bound is calculated from this chord, found with binary.
lval <- lowerbound(sampler,coefficients,index)

#The value of the upper bound is calculated by Wilson.
uval <- du.unnormalized(sampler, abscissae.result, z, lb, ub)
uniform <- runif(1)  


while(length(finalValues) < n){
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
      if(uniform <= exp(hValue - uval)){
        finalValues <- c(finalValues, sampler)
        abscissae.result = update(abscissae.result,sampler,hValue,h.deriv)
        coefficients = update_coeff(coefficients,sampler,abscissae.result)
        ############# Renee's function of check.log.concave #############
        if(check.log.concave(abscissae.result) == FALSE){
          stop("f is not log-concave!")
        }
        #################################################################
      }
    }
  }
}