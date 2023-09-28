#This file contains the execution code for the FSRMCD and IRMCD methods

#Function to execute the methods
calc_mcd <- function(z, method = "fsrmcd"){
  
  dados <- z$dadosc
  
  if(method == "fsrmcd"){
    
    tini <- Sys.time()
    
    outliers <- cerioli2010.fsrmcd.test(dados)$outliers
    
    rt <- as.numeric(difftime(Sys.time(), tini, units = "secs"))
    
  }else{
    
    tini <- Sys.time()
    
    outliers <- cerioli2010.irmcd.test(dados)$outliers
    
    rt <- as.numeric(difftime(Sys.time(), tini, units = "secs"))
    
  }
  
  #results table
  tab <- table(z$outver, outliers)
  
  if(nrow(tab) == 2){
    if(ncol(tab) == 2){
      
      ac <- tab %>% diag %>% sum %>% divide_by(nrow(dados))
      se <- tab[2,2] / tab %>% rowSums() %>% .[2]
      es <- tab[1,1] / tab %>% rowSums() %>% .[1] 
      
    }else{
      
      ac <- tab %>% diag %>% sum %>% divide_by(nrow(dados))
      se <- 0
      es <- tab[1,1] / tab %>% rowSums() %>% .[1]
      
    }
    
  }else{
    
    ac <- tab[1,1] %>% divide_by(nrow(dados))
    se <- NaN
    es <- tab[1,1] / tab %>% rowSums() %>% .[1]
    
  }
  
  return(data.frame(Sensivity = se, Specificity = es, Accuracy = ac, Runtime = rt, 
                    SQD = 0, SQE = 0, Ni_min = 0, K = 0, delta_est = 0, BIC = 0))
}