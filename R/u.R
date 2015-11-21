#### PUT ALL THE INPUT VARIABLES HERE ####

# f <- dnorm()
# h <- NULL
# k <- 100
# n <- 1000
#-----------------------------------------------------------------------------------------
# TEMPORARY INPUT VALUE
h <- function(x){
  return(log(dnorm(x)))
}
k <- 100
n <- 1000

abscissae.grid <- seq(-5, 5, length.out = 10)
gen.abscissae <- function(abscissae.grid, h){
  library(numDeriv)
  h.x <- h(abscissae.grid)
  h.deriv <- grad(func = h, x = abscissae.grid)
  abscissae.result <- cbind(abscissae.grid, h.x, h.deriv)
  return(abscissae.result)
}

abscissae.result <- gen.abscissae(abscissae.grid, h)
diff.h <- diff(abscissae.result[,2])
xh.deriv <- abscissae.result[,1] * abscissae.result[,3]
diff.xh.deriv <- diff(xh.deriv)
diff.h.deriv <- diff(abscissae.result[,3])

z <- (diff.h - diff.xh.deriv)/ (-diff.h.deriv)
#-----------------------------------------------------------------------------------------
all.mass <- sum(
  rep(c(-1, 1), length(abscissae.result[,3])) *
    rep(1/abscissae.result[,3], each = 2) * 
    exp(rep(abscissae.result[,2], each = 2) + 
          rep(abscissae.result[,3], each = 2) * 
          (c(-Inf, rep(z, each = 2), Inf) - rep(abscissae.result[,1], each = 2)))
)

norm.constant <- -log(all.mass)


du.normalized <- function(x, abscissae.result, z, norm.constant){
  j <- sapply(x, function(it) min(which(it <= c(z, Inf))))
  dens.u <- abscissae.result[j, 2] + (x - abscissae.result[j, 1]) * abscissae.result[j, 3] + norm.constant
  return(dens.u)
}


#-----------------------------------------------------------------------------------------

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

# cdf <- c(0,cum.mass)
# cdf <- seq(0,1, by = 0.1)
# x.hat <- S_inv(cdf, abscissae.result, z, norm.constant)


n.sim <- 1000
rs <- function(n.sim, S_inv, abscissae.result, z, norm.constant){
  u.temp <- runif(n = n.sim)
  x.temp <- S_inv(u.temp, abscissae.result, z, norm.constant)
  return(x.temp)
}

#-----------------------------------------------------------------------------------------
rs(n.sim, S_inv, abscissae.result, z, norm.constant)
x <- seq(-6, 6, by = 0.1)
du.normalized(x, abscissae.result, z, norm.constant)

haha <- function(x){ exp(du.normalized(x, abscissae.result, z, norm.constant))}
curve(haha, from = -3, to = 3)

lines(density(rs(1e6, S_inv, abscissae.result, z, norm.constant)), lty = 2 , col = 'red')

