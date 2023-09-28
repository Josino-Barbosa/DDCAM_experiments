#This file contains the code for the comparison experiments between the multivariate outlier 
#detection methods DDCAM, FSRMCD, IRMCD and MVE. These experiments are presented in the 
#article "Data-driven cluster analysis method: a novel outliers detection method in multivariate data"  

#Methods are tested on multivariate normal populations contaminated with outliers, 
#on multivariate exponential contaminated populations, and 
#on multivariate normal populations contaminated using point mass contamination

##Note: Since the authors are brazillian, some terms may be in Portuguese

#Packages and functions
library(pacman)
p_load("robustbase", "mvtnorm", "proxy", "MASS", "gencor", "magrittr", "ggplot2", "data.table", 
       tclust, gencor, googledrive, tidyverse, CerioliOutlierDetection)

#Sets the work directory
setwd(choose.dir())

#The codes for implementing the methods and generating contaminated data are in separate 
#files and will be loaded next

#Loads the function that generates contaminated multivariate normal populations
source("contaminate.R")

#Loads the function that generates contaminated multivariate exponential populations
source("contaminate_exp.R")

#Loads the function that generates contaminated multivariate normal populations using point mass contamination
source("contaminate_PMC.R")

#Loads the function that executes the DDCAM method
source("ddcam_max_kmeans.R")

#Loads the function that executes the MVE method
source("mve.R")

#Loads the function that executes the FSRMCD and IRMCD methods
source("metodos_mcd.R")


##Experiment 1 - Comparison of DDCAM, FSRMCD, IRMCD and MVE - multivariate normal data

#Random seed
set.seed(1234)

#Definition of simulation parameters
dim <- c(10, 20, 30, 50, 90)
delta <- c(0, 0.01, 0.025, 0.04, 0.05)
tam <- c(100, 200, 300, 400, 500)
nsim <- 200
cont <- 0

#Define objects to store the results
ddcam_normal <- NULL
mve_normal <- NULL
fsrmcd_normal <- NULL
irmcd_normal <- NULL

#Creates the seed vector of the replicas
vseed <- sample(1:(length(dim)*length(delta)*length(tam)*nsim), length(dim)*length(delta)*length(tam)*nsim,  replace = F)
v <- 1

for(j in delta){
  for(i in dim){
    for(h in tam){
      
      ini <- Sys.time()
      
      #Generate the correlation matrix
      Sig <- gencor(i, method = "medium", lim_low = 0.2, lim_medium = 0.8)$Matrix
      
      for(k in 1:nsim){
        
        print(paste(i, j, k, h, sep = " - "))
        
        #Fix the seed of replicas
        set.seed(vseed[v])
        
        #Data generation
        z <- contaminate(n = h, delta = j, mu1 = rnorm(i, 0, 1), mu2 = rnorm(i, 2, 1), sigma = Sig)
        
        #Application of methods
        
        #DDCAM
        ddcam_normal <- rbind(ddcam_normal, 
                                     c("ddcam", i, j, h, ddcam_max(z, sd_dist = 2.0)))
        
        #FSRMCD
        fsrmcd_normal <- rbind(fsrmcd_normal, 
                                     c("fsrmcd", i, j, h, calc_mcd(z, method = "fsrmcd")))
        
        #IRMCD
        irmcd_normal <- rbind(irmcd_normal, 
                               c("irmcd", i, j, h, calc_mcd(z, method = "irmcd")))
        
        #MVE
        mve_normal <- rbind(mve_normal, 
                              c("mve", i, j, h, mve(z)))
        
        #Increment for changing the seed of replicas
        v <- v+1
        
      }
    }
  }
}

#Experiment 1 results
resuluts_normal <- rbind(ddcam_normal, fsrmcd_normal, irmcd_normal, mve_normal) %>% 
  data.frame %>% 
  rename("metodo" = "V1", "dim" = V2, "delta" = V3, "n" = V4) %>% 
  arrange(dim, delta, n, metodo) %>% 
  mutate(metodo = metodo %>% as.character,
         dim = dim %>% as.numeric,
         delta = delta %>% as.numeric,
         n = n %>% as.numeric,
         Sensivity = Sensivity %>% as.numeric,
         Specificity = Specificity %>% as.numeric,
         Accuracy = Accuracy %>% as.numeric,
         Runtime = Runtime %>% as.numeric,
         SQD = SQD %>% as.numeric,
         SQE = SQE %>% as.numeric,
         Ni_min = Ni_min %>% as.numeric,
         K = K %>% as.numeric,
         delta_est = delta_est %>% as.numeric,
         BIC = BIC %>% as.numeric) %>% 
  group_by(delta, dim, n, metodo) %>% 
  summarise(ac = mean(Accuracy), 
            se = mean(Sensivity), 
            sp = mean(Specificity), 
            rt = mean(Runtime))


##Experiment 2 - Comparison DDCAM, FSRMCD, IRMCD and MVE - multivariate exponential data

#Random seed
set.seed(1234)

#Definition of simulation parameters
dim <- c(10, 20, 30, 50, 90)
delta <- c(0, 0.01, 0.025, 0.04, 0.05)
tam <- c(100, 200, 300, 400, 500)
nsim <- 200
cont <- 0

#Define objects to store the results
ddcam_exp <- NULL
mve_exp <- NULL
fsrmcd_exp <- NULL
irmcd_exp <- NULL

#Creates the seed vector of the replicas
vseed <- sample(1:(length(dim)*length(delta)*length(tam)*nsim), length(dim)*length(delta)*length(tam)*nsim,  replace = F)
v <- 1

for(j in delta){
  for(i in dim){
    for(h in tam){
      
      ini <- Sys.time()
      
      #Generate the correlation matrix
      Sig <- gencor(i, method = "medium", lim_low = 0.2, lim_medium = 0.8)$Matrix
      
      for(k in 1:nsim){
        
        print(paste(i, j, k, h, sep = " - "))
        
        #Fix the seed of replicas
        set.seed(vseed[v])
        
        #Data generation
        z <- contaminate_exp(nsample = h, delta = j, mu1 = rep(1, i), mu2 = 2, sigma2 = Sig)
        
        #Application of methods
        
        #DDCAM
        ddcam_exp <- rbind(ddcam_exp, 
                              c("ddcam", i, j, h, ddcam_max(z, sd_dist = 2.0)))
        
        #FSRMCD
        fsrmcd_exp <- rbind(fsrmcd_exp, 
                               c("fsrmcd", i, j, h, calc_mcd(z, method = "fsrmcd")))
        
        #IRMCD
        irmcd_exp <- rbind(irmcd_exp, 
                              c("irmcd", i, j, h, calc_mcd(z, method = "irmcd")))
        
        #MVE
        mve_exp <- rbind(mve_exp, 
                            c("mve", i, j, h, mve(z)))
        
        #Increment for changing the seed of replicas
        v <- v+1
        
      }
    }
  }
}

##Experiment 2 results
resuluts_exp <- rbind(ddcam_exp, fsrmcd_exp, irmcd_exp, mve_exp) %>% 
  data.frame %>% 
  rename("metodo" = "V1", "dim" = V2, "delta" = V3, "n" = V4) %>% 
  arrange(dim, delta, n, metodo) %>% 
  mutate(metodo = metodo %>% as.character,
         dim = dim %>% as.numeric,
         delta = delta %>% as.numeric,
         n = n %>% as.numeric,
         Sensivity = Sensivity %>% as.numeric,
         Specificity = Specificity %>% as.numeric,
         Accuracy = Accuracy %>% as.numeric,
         Runtime = Runtime %>% as.numeric,
         SQD = SQD %>% as.numeric,
         SQE = SQE %>% as.numeric,
         Ni_min = Ni_min %>% as.numeric,
         K = K %>% as.numeric,
         delta_est = delta_est %>% as.numeric,
         BIC = BIC %>% as.numeric) %>% 
  group_by(delta, dim, n, metodo) %>% 
  summarise(ac = mean(Accuracy), 
            se = mean(Sensivity), 
            sp = mean(Specificity), 
            rt = mean(Runtime))


##Experiment 3 - Comparison of DDCAM, FSRMCD, IRMCD and MVE - Normal data using point mass contamination

#Random seed
set.seed(1235)

#Definition of simulation parameters
dim <- c(10, 20, 30)
delta <- c(0.025, 0.05, 0.075, 0.10, 0.15)
tam <- c(100, 200, 300, 400, 500)
nsim <- 200
cont <- 0

#Define objects to store the results
ddcam_PMC <- NULL
mve_PMC <- NULL
fsrmcd_PMC <- NULL
irmcd_PMC <- NULL

#Creates the seed vector of the replicas
vseed <- sample(1:(length(dim)*length(delta)*length(tam)*nsim), length(dim)*length(delta)*length(tam)*nsim,  replace = F)
v <- 1

for(j in delta){
  for(i in dim){
    for(h in tam){
      
      ini <- Sys.time()
      
      #Generate the correlation matrix
      Sig <- gencor(i, method = "medium", lim_low = 0.2, lim_medium = 0.8)$Matrix
      
      for(k in 1:nsim){
        
        print(paste(i, j, k, h, sep = " - "))
        
        #Fix the seed of replicas
        set.seed(vseed[v])
        
        #Data generation
        z <- contaminate_PMC(n = h, delta = j, mu1 = rnorm(i, 0, 1), sigma = Sig)

        #Application of methods
        
        #DDCAM
        ddcam_PMC <- rbind(ddcam_PMC, 
                              c("ddcam", i, j, h, ddcam_max(z, sd_dist = 2.0)))
        
        #FSRMCD
        fsrmcd_PMC <- rbind(fsrmcd_PMC, 
                               c("fsrmcd", i, j, h, calc_mcd(z, method = "fsrmcd")))
        
        #IRMCD
        irmcd_PMC <- rbind(irmcd_PMC, 
                              c("irmcd", i, j, h, calc_mcd(z, method = "irmcd")))
        
        #MVE
        mve_PMC <- rbind(mve_PMC, 
                            c("mve", i, j, h, mve(z)))
        
        #Increment for changing the seed of replicas
        v <- v+1
        
      }
    }
  }
}

##Experiment 3 results
resuluts_pmc <- rbind(ddcam_PMC, fsrmcd_PMC, irmcd_PMC, mve_PMC) %>% 
  data.frame %>% 
  rename("metodo" = "V1", "dim" = V2, "delta" = V3, "n" = V4) %>% 
  arrange(dim, delta, n, metodo) %>% 
  mutate(metodo = metodo %>% as.character,
         dim = dim %>% as.numeric,
         delta = delta %>% as.numeric,
         n = n %>% as.numeric,
         Sensivity = Sensivity %>% as.numeric,
         Specificity = Specificity %>% as.numeric,
         Accuracy = Accuracy %>% as.numeric,
         Runtime = Runtime %>% as.numeric,
         SQD = SQD %>% as.numeric,
         SQE = SQE %>% as.numeric,
         Ni_min = Ni_min %>% as.numeric,
         K = K %>% as.numeric,
         delta_est = delta_est %>% as.numeric,
         BIC = BIC %>% as.numeric) %>% 
  group_by(delta, dim, n, metodo) %>% 
  summarise(ac = mean(Accuracy), 
            se = mean(Sensivity), 
            sp = mean(Specificity), 
            rt = mean(Runtime))