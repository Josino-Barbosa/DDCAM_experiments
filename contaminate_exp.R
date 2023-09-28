#This file generates contaminated multivariate exponential populations

#Install the lcmix package if necessary
#install.packages("lcmix", repos=c("http://R-Forge.R-project.org",
#                                  "http://cran.at.r-project.org"),dependencies=TRUE)

#Used packages
library(lcmix)

#Function to contaminate a multivariate exponential population
contaminate_exp <- function(nsample, delta, mu1, mu2, sigma2, sort = T){
  
  #For contaminations that would generate less than 1 outlier, the contamination
  #rate becomes the probability of generating an outlier drawn via Bernoulli
  
  fat_cont <- nsample*delta
  
  if(fat_cont < 1 & delta != 0){
    
    delta = rbinom(1, 1, fat_cont)/nsample
    
  }
  
  if(delta != 0){
    if(round(nsample*delta) == 1){
      
      dadosc <- rbind(rmvexp(n = 2*round(nsample*delta), rate = mu1, corr = sigma2)[1, ],
                      rmvexp(n = nsample - round(nsample*delta), rate = mu1, corr = sigma2) + 
                        matrix(mu2, nrow = nsample - round(nsample*delta), ncol = ncol(sigma2)))
      
    }else{
      
      dadosc <- rbind(rmvexp(n = round(nsample*delta), rate = mu1, corr = sigma2) + 
                        matrix(mu2, nrow = round(nsample*delta), ncol = ncol(sigma2)),
                      rmvexp(n = nsample - round(nsample*delta), rate = mu1, corr = sigma2))
      
    }
    
    outver <- c(rep(1, round(nsample*delta)), rep(0, nsample - round(nsample*delta)))
  }
  
  if(delta == 0){
    dadosc <- rmvexp(nsample, rate = mu1, corr = sigma2)
    outver <- rep(0, nsample)
  }
  
  if(sort == T){
    id <- sample(1:nsample, nsample, replace = F)
    dadosc <- dadosc[id, ]
    outver <- outver[id]
  }
  
  invisible(list(dadosc = dadosc, outver = outver))
}
