#---------------------------------------------
#Set up arguments:
n <- 1000

#IMPORTANT!  Note on Tk:  I picked 0 and two values just to the left and right of 0,
#then updated my own Tk by assining any test values of x to become new Tk values 
#if they fall outside the boundries.  Then I draw new chords.  Obviously the 'updating'
#part of the package will replace this updating method, it was just so I could run
#my code.

Tk <- c(-0.001, 0, 0.001)
#---------------------------------------------

#Example function for f
h <- function(x){
  return(log(dt(x, df = 10)))}

#---------------------------------------------
#This is the main engine.  It takes the log of a function (h), a required number of
#values (n), and the starting points to draw chords (the vector Tk).  
level2 <- function(h, n, Tk){
  #'finalValues' will be what we return.  We initiate it empty.
  finalValues <- c()
  
  #'Coefficients' are the sope and y intercept associated with each chord in the 'l'
  #function  They are in a matrix that gets created by 'lupdater'.
  coefficients <- lupdater(Tk)
  
  #Run the code until we hit the desired length
  while(length(finalValues) < n){
    #Take a random point from the 'sk' function, known as 'rs'.
    sampler <- (rs(1, S_inv, abscissae.result, z, norm.constant))
    #If that point lies outsidethe bounds set by my Tk values
    if(sampler < Tk[1] | sampler > Tk[length(Tk)]){
      #Add that point to become a new Tk value.
      Tk <- sort(c(sampler, Tk))
      #And update the coefficient matrix (since we added a new chord)
      coefficients <- lupdater(Tk)
    
    }else{ 
      #First we run a binary search to find which chord the point finds itself within.
      index <- binary(Tk, sampler)
      
      #The value of the lower bound is calculated from this chord, found with binary.
      lval <- lowerbound(sampler,coefficients,index)
      
      #The value of the upper bound is calculated by Wilson.
      uval <- du.normalized(sampler, abscissae.result, z, norm.constant)
      
      #We also grab a number from the uniform distribution.
      uniform <- runif(1)  
      
      #The big IF statements, this one's the 'squeezing' step
      if(uniform <= exp(lval - uval)){
        finalValues <- c(finalValues, sampler)
        
      }else{
        #Only evaluate the log of the function if we fail to squeeze.  Rejection step.
        hValue <- h(sampler)
        if(uniform <= exp(hValue - uval)){
          finalValues <- c(finalValues, sampler)
          
          #Also update this value into Tk.  Sorry if this gets done in a later step, 
          #we can always remove this functionality!
          Tk <- sort(c(sampler, Tk))
          coefficients <- lupdater(Tk)
          
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

#This function updates the chords of the 'l' function with new Tk points.
lupdater <- function(Tk){
  coefficients <- matrix(, nrow = (length(Tk)-1), ncol = 2)
  for(j in 1:(length(Tk)-1)){
    coefficients[j,1] <- (h(Tk[j])-h(Tk[j+1]))/(Tk[j]-Tk[j+1])   
    coefficients[j,2] <- h(Tk[j]) - coefficients[j,1]*Tk[j]
  }
  return(coefficients)
}
    


