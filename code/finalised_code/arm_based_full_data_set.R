# Packages
pacman::p_load(metafor, MCMCglmm, tidyverse, rotl, phytools, corrplot, ape)
library(dplyr)
library(tidyr)
library(readr)
library(purr)
library(tibble)
library(stringr)
library(forcats)
library(ggplot2)


source("./code/func.R")

#Data processing
rerun = FALSE

if(rerun == TRUE){
  
  # Bring in the data and convert key variables to required classes. 
  data <- read.csv("./data/finalised_data_sets/arm_based_full_dataset.csv", stringsAsFactors = FALSE)
  data$body_mass_g <- as.numeric(data$body_mass_g)
  data$study_ID <- as.factor(data$study_ID)
  data$study_type <- as.factor(data$study_type)
  data$respiration_mode <- as.factor(data$respiration_mode)
 
  
  # Separating genus and species
  genus <- as.data.frame(do.call("rbind", str_split(str_trim(data$species, side = "both"), " ")))
  names(genus) <- c("genus", "species_new")
  data <- cbind(data, genus)
  data$species_rotl <- paste0(data$genus, "_", data$species_new)
  data <- data[,-22]	
  
  
  # Tcentering 
  #Tw-centering of tempreature within each species across all studies on that species
  sp <- split(data, data$species)
  data$T_w <- do.call("rbind", lapply(sp, function(x) scale(x$T, scale = FALSE)))
  
  # Add in observation level random effects
  data$obs <- 1:dim(data)[1]
  
  # Write the file, so it can be loaded more easily 
  write.csv(data, "./data/data_arm.csv", row.names=FALSE)
  
}else {
  
  data <- read.csv("./data/data_arm.csv")
  
}


#Importing phylogeny from TimeTree
  # Import TimeTree phylogeny
  tree <- read.tree("./data/full_phylogeny.NWK")
  plot(tree)
  A <- inverseA(tree, nodes = "TIPS")$Ainv


  # Check what is different-- different number of species (16 compared to 13)
  setdiff(unique(data$species_rotl), sort(rownames(A)))

  # Check same number of levels
  length(rownames(A))
  length(unique(data$species_rotl))
  
  
# Calculate the sampling variance for all mean estimates
  data$mean_sv <- with(data, m_sv(mean, sd, n))
  data$sd_sv <- with(data, sd_sv(n))
  
 
  
# Bayesian Priors              
  data2 <- data[complete.cases(data[,c("T_w", "body_mass_g", "respiration_mode")]),]
  
  prior_slope <- list(R = list(V = 1, nu = 0.001),
                      G = list(G1 = list(V=1, nu = 0.02),
                               G2 = list(V = diag(2), nu = 2)))  
  
  
  
     # Exploratory analysis 
  with(data, plot(log(mean) ~ log(sd)))
  with(data, hist(log(mean)))
  with(data, hist(log(sd)))
  ggplot(data, aes(y = log(mean), x=T_w, colour= species))+
    geom_line()
  
  data %>% filter(mean >100)
  data %>% filter(species_rotl == "Crocodylus_porosus")
  
  # Have a look at some details about the temperatures
  Temps <- data2 %>% group_by (species_rotl) %>% summarise(Num_T = length(unique(T)), rangeL = range(T)[1], rangeU = range(T)[2], diff = range(T)[2] - range(T)[1])
  range(Temps$diff)

  
  
  ###########################################
  #Mean models 
  ###########################################	
  
  model1.mean <- MCMCglmm(
    log(mean) ~  T_w + log(body_mass_g) + respiration_mode * T_w, 
    mev = data2$mean_sv, 
    random = ~us(1):study_ID + us(1 + T_w):species_rotl, 
    ginverse = list(species_rotl = A), 
    data = data2, 
    prior = prior_slope, 
    nitt = 130000, burnin = 30000, thin = 50)
  summary(model1.mean)
  plot(model1.mean)
  
  # Get the posterior estimate for temperature effect within species
  T_w <- model1.mean$Sol[,"T_w"]
  mean(T_w)
  HPDinterval(T_w)
  ((1-(exp(mean(T_w))))*7.4)*100
  
  
  ###########################################
  #Variance models – CVR models
  ###########################################
  
  #SD model
  model2.sd <- MCMCglmm(
    log(sd) ~ log(mean) + T_w + log(body_mass_g) + respiration_mode,
    mev = data2$sd_sv, 
    random = ~us(1):study_ID + us(1 + T_w):species_rotl, 
    ginverse = list(species_rotl = A), 
    data = data2, 
    prior = prior_slope, 
    nitt = 130000, burnin = 30000, thin = 50)
  summary(model2.sd)
  plot(model2.sd)
  
  # Get the posterior estimate for temperature effect within species
  T_w_sd <- model2.sd$Sol[,"T_w"]
  mean(T_w_sd)
  HPDinterval(T_w_sd)
  
