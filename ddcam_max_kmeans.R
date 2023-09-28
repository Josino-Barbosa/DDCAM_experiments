#This file contains the execution code for the DDCAM method

##Note: Since the authors are brazillian, some terms may be in Portuguese

#Used packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(proxy, robustbase, magrittr)

#DDCAM function
ddcam_max <- function(dadosc, sd_dist = 2.0, tau = 2, delta = NULL){

  dadosc$dadosc <- scale(dadosc$dadosc)
  
  #Delta calculation
  sddist <- function(dados){
    sdcol <- tau*sd(dados)
    cont_dist <- function(x, sdcol){
      cd <- sum(abs(x - mean(x)) > sdcol)
      return(cd)
    }
    
    k <- sum(cont_dist(dados, sdcol))
    return(k)
    
  }
  
  delta_est <- max(apply(dadosc$dadosc, 2, sddist)/nrow(dadosc$dadosc))
  
  if(is.null(delta) == T){
    delta <- delta_est
  }
  
  #Start runtime
  t1 <- Sys.time()
  
  resultados_gerais <- NULL
  
  for(nclust in 2:(ceiling(nrow(dadosc$dadosc)/log(nrow(dadosc$dadosc))) + 1)){
    
    #Number of centroids
    k <- nclust
    
    ###Choose the initial centroid for the k-means using the Ward method
    
    #Generates the grouping
    hc <- hclust(dist(dadosc$dadosc), method = 'ward.D2')
    
    #Identifies individuals by group
    id_grup <- cutree(hc, k)
    
    #Calculate the centroids
    cent <- NULL
    
    for(i in 1:k){
      
      cent <- rbind(cent, colMeans(dadosc$dadosc[id_grup == i, , drop = FALSE]))
      
    }
    
    #Aggregates data and outlier identifier
    pop <- cbind(dadosc$dadosc, dadosc$outver)
    
    #Generate K-means clustering
    KM <- kmeans(dadosc$dadosc, centers = cent)
    
    #Counter of at least one success
    contref <- 0
    
    #First Refinement
    if(sum(KM$size <= delta * sum(KM$size)) != 0){
      
      #Calculation of centroids
      dk <- KM$centers
      
      #Identification of individuals by group
      cl <- KM$cluster
      dadosA <- cbind(pop, cl)
      
      #Data median
      medianA <- colMedians(dadosc$dadosc)
      
      #Distance between the centroids and the median
      d <- dist(rbind(dk, medianA))
      d <- as.matrix(d)
      di <- d[(k+1), -(k+1)]
      
      #Defines whether there are groups of outliers
      id_clust <- as.numeric(di > sd_dist * sd(di))
      
      #Checks if groups are smaller than n*delta
      id_size <- as.numeric(KM$size <= delta * sum(KM$size))
      
      #Second Refining
      if(sum(id_clust * id_size) > 0){
        
        #Counter of at least one success
        contref <- 1
        out_detect <- as.numeric(cl %in% which((id_clust * id_size) == 1))
        
        results <- data.frame(Detected = factor(out_detect, levels = c(0,1), labels = c(0,1)), 
                              Outlier = factor(dadosc$outver, levels = c(0,1), labels = c(0,1)))
        tab <- table(results)
        tab <- tab[,2:1]
        tab <- tab[2:1,]
        
        ###########################################################
        
        resultados_gerais <- rbind(resultados_gerais, data.frame(Sensivity = tab[1, 1] / (tab[1, 1] + tab[2, 1]),
                                                                 Specificity = tab[2, 2] / (tab[1, 2] + tab[2, 2]),
                                                                 Accuracy = (tab[1, 1] + tab[2, 2]) / sum(tab),
                                                                 Runtime = NA,
                                                                 SQD = KM$tot.withinss,
                                                                 SQE = KM$betweenss,
                                                                 Ni_min = min(KM$size),
                                                                 K = k,
                                                                 delta_est = delta))
        ###########################################################
        
      } #End of second refinement
    }#End of the first refinement
  
  }
  
  #Closes execution time
  t2 <- Sys.time()
  
  #Calculates execution time
  runtime <- as.numeric(difftime(t2, t1, units = "secs"))
  
  resultados_gerais$Runtime <- runtime
  
  #Calculation of performance measures for zero detections
  if(contref == 0){
    
    out_detect <- rep(0, nrow(dadosc$dadosc))
    
    #Closes execution time
    t2 <- Sys.time()
    
    results <- data.frame(Detected = factor(out_detect, levels = c(0,1), labels = c(0,1)), 
                          Outlier = factor(dadosc$outver, levels = c(0,1), labels = c(0,1)))
    tab <- table(results)
    tab <- tab[,2:1]
    tab <- tab[2:1,]
    
    ###########################################################
    #Calculates execution time
    runtime <- as.numeric(difftime(t2, t1, units = "secs"))
    
    #Results when outliers are not detected
    resultados_gerais <- rbind(resultados_gerais, data.frame(Sensivity = tab[1, 1] / (tab[1, 1] + tab[2, 1]),
                                                             Specificity = tab[2, 2] / (tab[1, 2] + tab[2, 2]),
                                                             Accuracy = (tab[1, 1] + tab[2, 2]) / sum(tab),
                                                             Runtime = runtime,
                                                             K = NA,
                                                             delta_est = delta))
    ###########################################################
    
  }
  
  #Calculation of the best K via BIC
  n <- nrow(dadosc$dadosc)
  d <- ncol(dadosc$dadosc)
  K <- resultados_gerais$K
  sqd <- resultados_gerais$SQD %>% as.numeric
  resultados_gerais$BIC <- sqd + log(n)*K*d
  
  #Results if outliers are detected
  if(is.na(resultados_gerais$BIC)[1] != T){
    return(resultados_gerais[which.min(resultados_gerais$BIC), ])
  }
  
  #Results if no outliers are detected
  if(is.na(resultados_gerais$BIC)[1] == T){
    return(resultados_gerais %>% data.frame)
  }
  
}