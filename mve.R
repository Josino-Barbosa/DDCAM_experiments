#This file contains the execution code for the MVE method

mve <- function(dadosc, conf = 0.05){
  
  #Start runtime
  t1 <- Sys.time()
  
  #Number of rows and columns in the data
  n <- nrow(dadosc$dadosc)
  p <- ncol(dadosc$dadosc)
  
  #MEV estimates of location and scale parameters
  estrob <- cov.mve(dadosc$dadosc, 
                    cor = T,
                    method = "mve")
  
  #Rental Parameter
  Tx <- estrob$center
  
  #Scale parameter
  Cx <- estrob$cov
  
  #Robust mahalanobis distance via MVE
  invCx <- solve(Cx)
  RD <- matrix(0, n, 1)
  quadrRD <- RD
  
  for(i in 1:n){
    RD[i] <- ((t(dadosc$dadosc[i, ] - Tx)) %*% invCx %*% (dadosc$dadosc[i, ] - Tx)) ^ .5
    quadrRD[i] <- (RD[i])^2
  }
  
  #Chi-square quantile to detect multivariate outliers
  confianca <- 1 - conf
  quantil <- qchisq(confianca, p)
  
  #Creates the matrix that indicates whether the observation is an outlier
  out_detect <- matrix(0, n, 1)
  
  for(i in 1:n){
    if(quadrRD[i] > quantil){
      out_detect[i] <- 1
    }
    else{
      out_detect[i] <- 0
    }
  }
  
  #Closes execution time
  t2 <- Sys.time()
  
  results <- data.frame(Detected = factor(out_detect, levels = c(0,1), labels = c(0,1)), 
                        Outlier = factor(dadosc$outver, levels = c(0,1), labels = c(0,1)))
  tab <- table(results)
  tab <- tab[,2:1]
  tab <- tab[2:1,]
  
  #Calculates execution time
  runtime <- as.numeric(difftime(t2, t1, units = "secs"))
  
  return(data.frame(Sensivity = tab[1, 1] / (tab[1, 1] + tab[2, 1]),
    Specificity = tab[2, 2] / (tab[1, 2] + tab[2, 2]),
    Accuracy = (tab[1, 1] + tab[2, 2]) / sum(tab),
    Runtime = runtime, 
    SQD = 0, SQE = 0, Ni_min = 0, K = 0, delta_est = 0, BIC = 0))
}
