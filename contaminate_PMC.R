#This file generates contaminated multivariate normal populations using point mass contamination

#Used packages
library(mvtnorm, dplyr)

#Function to contaminate a multivariate normal population
contaminate_PMC <- function(n, delta, mu1, sigma2, sort = T){
  
  #For contaminations that would generate less than 1 outlier, the contamination
  #rate becomes the probability of generating an outlier drawn via Bernoulli
  
  fat_cont <- ceiling(n*delta)
  
  if(delta != 0){
    dados <- rmvnorm(n, mean = mu1, cov = sigma2)
    
    cont_row <- (nrow(dados) - fat_cont + 1 ):nrow(dados)
    
    mui <- mu1
    sdi <- diag(cov(dados))
    
    point_contamination <- mui + sample(c(-2, 2), ncol(dados), replace = T) * sdi
    
    contaminated_rows <- rep(point_contamination, fat_cont) %>% 
      matrix(nrow = fat_cont, byrow = T) 
    
    dadosc <- rbind(dados[-cont_row, ], contaminated_rows)
    
    outver <- numeric(nrow(dados))
    outver[cont_row] <- 1
    
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