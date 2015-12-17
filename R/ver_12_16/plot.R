
set.seed(0)
k=30
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
lvalf <- function(sampler) lowerbound(sampler,coefficients,index)

#The value of the upper bound is calculated by Wilson.
uvalf <- function(sampler) du.unnormalized(sampler, abscissae.result, z, lb, ub)

first_step <- function(x) exp(lvalf(x) - uvalf(x))
second_step <- function(x) exp(h(x) - uvalf(x))

library(ggplot2)

pl = ggplot(data.frame(x=c(0,1)), aes(x)) +
  stat_function(fun=first_step, geom="line", aes(colour="phase 1")) +
  stat_function(fun=second_step, geom="line", aes(colour="phase 2")) 
 


 


x = c()
y = c()
Z = c()
for( i in 1:10){
  z = compute_z(abscissae.result)
  norm.constant = compute_norm_constant(abscissae.result,z, lb, ub)
  sampler <- rs(1, S_inv, abscissae.result, z, norm.constant, lb, ub)[1,1]
  uniform <- runif(1)  
  if(sampler < Tk[1] | sampler > Tk[length(Tk)]){
    
    z = "rejected"
  
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
      z = "phase 1"
    }else{
      # Only evaluate the log of the function if we fail to squeeze.  Rejection step.
      hValue <- h(sampler)
      h.deriv <- grad(func = h, x = sampler) 
      if(uniform <= exp(hValue - uval)){ 
        z = "phase 2" 
        }else{
        z = "rejected"
        }
        #################################################################
    }
  }
  x = c(x,sampler)
  y = c(y,uniform)
  Z = c(Z,z)
  
}
 
df = data.frame(x, y, Z) 

pl + geom_point(data=df, aes(x = x, y = y,colour = Z))  











