#This file generates contaminated multivariate normal populations

#Used packages
library(mvtnorm)

#Function to contaminate a multivariate normal population
contaminate <- function(n, delta, mu1, mu2, sigma2, sort = T){
  
  #For contaminations that would generate less than 1 outlier, the contamination
  #rate becomes the probability of generating an outlier drawn via Bernoulli
  
  fat_cont <- n*delta
  
  if(fat_cont < 1 & delta != 0){
    
    delta = rbinom(1, 1, fat_cont)/n
    
  }
  
  if(delta != 0){
    dadosc <- rbind(rmvnorm(round(n*delta), mean = mu1, cov = sigma2),
                    rmvnorm(n - round(n*delta), mean = mu2, cov = sigma2))
    outver <- c(rep(1, round(n*delta)), rep(0, n - round(n*delta)))
  }
  
  if(delta == 0){
    dadosc <- rmvnorm(n, mean = mu1, cov = sigma2)
    outver <- rep(0, n)
  }
  
  if(sort == T){
    id <- sample(1:n, n, replace = F)
    dadosc <- dadosc[id, ]
    outver <- outver[id]
  }
  
  invisible(list(dadosc = dadosc, outver = outver))
}