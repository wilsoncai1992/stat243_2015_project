n.sim <- 10^4
nu <- 10
#var <- 1

c <- gamma(11/2)*sqrt(pi)/sqrt(10)/gamma(5)*2*(1.1^(-11/2))

F.density <- function(x, nu = 10){
  result <- dt(x, df = nu) 
  return(result)
}

G <- function(n.sim){
  g.sample <- rcauchy(n = n.sim, location = 0, scale = 1)
  return(g.sample)
}
G.density <- function(x){
  result <- dcauchy(x = x, location = 0, scale = 1)
  return(result)
}


rejection.fn <- function(n.sim, f.density, g.simu, g.density, c ){
  ##### STEP 1 #####
  u.1 <- runif(n.sim)
  y.sample <- g.simu(n.sim)
  
  ##### STEP 2 #####
  x.output <- c()
  for(i in 1:n.sim){
    if(u.1[i] <= f.density(y.sample[i])/(c*g.density(y.sample[i])) ){
      x.output <- cbind(x.output,y.sample[i])
    }  
  }
  x.output <- as.vector(x.output)
  return(x.output)
}

####### TARGET SAMPLE SIZE = 10,000 #######
answer.vec <- rejection.fn(n.sim = n.sim*1.5, f = F.density, g.simu = G, g.density = G.density, c = c )
answer.vec <- tail(answer.vec, n = n.sim)
plot(density(answer.vec), lwd=2, main=expression(paste('Density Plot of True and Simulated Student-',t[10](0,1), ' Distribution')))

# true.vec <- rt(n = 10^5, df = nu)
# lines(density(true.vec), col='blue', lwd=2,lty=2)

true.density <- function(grid.sim = 10^2,nu = 10, var = 1){
  x.vec <- seq(-5,5,length.out = grid.sim)
  
  true.vec <- F.density(x = x.vec)
  
  return(list(true.vec, x.vec))
}

true.ans <- true.density()
lines(true.ans[[1]]~true.ans[[2]], col='blue',lty=2,lwd =2)

legend('topright', legend = c('True','Rejection Sampling'), col=c('blue','black'),lty=c(2,1),lwd=2)


# fun <- function(x){
#   f1 <- dt(x = x,df = 10)
#   g1 <- dcauchy(x = x,location = 0,scale = 1)
#   return(f1/g1)
# }
# 
# fun(seq(-5,5,length.out = 10^5))
# 
# max(fun(seq(-5,5,length.out = 10^6)))
# c
# gamma(11/2)*sqrt(pi)/gamma(5)/sqrt(10)
