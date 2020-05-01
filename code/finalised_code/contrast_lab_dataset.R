# Clean workspace
rm(list = ls()) 


# Load libraries
pacman::p_load(metafor, MCMCglmm, tidyverse, rotl, phytools, corrplot, devtools)
install_github("daniel1noble/metaAidR"); library(metaAidR)
source("./code/finalised_code/func.R")
library(dplyr)
library(tidyr)
library(readr)
library(purr)
library(tibble)
library(stringr)
library(forcats)

#data processing
rerun_data = FALSE
if(rerun_data == TRUE){
  
  # Bring in the data and convert key variables to required classes. 
  data <- read.csv("./data/finalised_data_sets/contrast_lab_dataset.csv", stringsAsFactors = FALSE)
  data$body_mass <- as.numeric(data$body_mass)
  data$t_magnitude <- ordered(data$t_magnitude, levels = c("plus3", "plus5-7", "plus8-9", "plus10"))
  data$study_ID <- as.factor(data$study_ID)
 
  
  # Separate genus and species
  genus <- as.data.frame(do.call("rbind", str_split(str_trim(data$species, side = "both"), " ")))
  names(genus) <- c("genus", "species_new", "sub_spp")
  data <- cbind(data, genus)
  
  # Fix up the species names so they match with phylogeny
  data$species_rotl <- paste0(data$genus, "_", data$species_new)
  
  
  # Write the file, so it can be loaded more easily 
  write.csv(data, "./data/data_full.csv", row.names=FALSE)
}else {
  data_lab <- read.csv("./data/data_full.csv", stringsAsFactors=FALSE)
  data_lab$t_magnitude <- ordered(data_lab$t_magnitude, 
                                    levels = c("plus3", "plus5-7", "plus8-9", "plus10"))
}




# Import TimeTree phylogeny
tree <- read.tree("./data/finalised_data_sets/phylo_lab.NWK")
plot(tree)



# Create a phylogenetic correlation matrix
PhyloA <- vcv(tree, corr = TRUE)



#Tree checks_full dataset
  setdiff(unique(data_lab$species_rotl), sort(rownames(PhyloA)))
   # Check same number of levels
  length(rownames(PhyloA))
  length(unique(data_lab$species_rotl))
 



# Have a look at the mean-variance relationship
plot_func(data_lab, "sd_t1", "mean_t1")

# Calculating effect sizes, measure="ROM" means lnRR es, append=true means will add es to data set. yi= es estimate, vi sampling variance.
data_lab_ROM <- escalc(m1i = mean_t2, 
                         m2i = mean_t1, 
                         n1i = n_t2, 
                         n2i = n_t1, 
                         sd1i = sd_t2, 
                         sd2i = sd_t1, 
                         append = TRUE, measure ="ROM", 
                         data = data_lab)

data_lab_CVR <- escalc(m1i = mean_t2, 
                         m2i = mean_t1, 
                         n1i = n_t2, 
                         n2i = n_t1, 
                         sd1i = sd_t2, 
                         sd2i = sd_t1, 
                         append = TRUE, measure ="CVR", 
                         data = data_lab)

# Error calc-residual variance at the observation level as metafor does not add this by default
data_lab_ROM$obs <- 1:dim(data_lab_ROM)[1]
data_lab_CVR$obs <- 1:dim(data_lab_CVR)[1]

# Make shared control matrix for ROM dataset
  data_lab_ROM$sc_cluster <- interaction(data_lab_ROM$study_ID, 
                                         data_lab_ROM$shared_control)

  # Create the Shared control V matrix
  V <- metaAidR::make_VCV_matrix(data_lab_ROM, "vi", "sc_cluster", 
                               type = "vcv", rho = 0.5)
  
  V <- make_VCV_matrix(data_lab_ROM, "vi", "sc_cluster", 
                                 type = "vcv", rho = 0.5)
  
  # Look at the shared control (co)variance matrix
  corrplot(as.matrix(V), is.corr = FALSE)
  write.csv(V, file = "./output_matrices/sc_matrix_rom.csv")


  # Make shared control matrix for CVR dataset
  data_lab_CVR$sc_cluster <- interaction(data_lab_CVR$study_ID, 
                                         data_lab_CVR$shared_control)
  # Create the SC V matrix
  V2 <- make_VCV_matrix(data_lab_CVR, "vi", "sc_cluster", type = "vcv", rho = 0.5)

  # Look at the shared control (co)variance matrix
  corrplot(as.matrix(V2), is.corr = FALSE)
  write.csv(V2, file = "./output_matrices/sc_matrix_cvr.csv")

  #Calculing weighted means and SD for delta_t
  weighted.mean(data_lab_ROM$ delta_t, data_lab_ROM$ n_t1)
  weightedSD(data_lab_ROM$delta_t, data_lab_ROM$n_t1)
  
  #Calculing weighted means and SD for body mass
  weighted.mean(data_lab_ROM$ body_mass, data_lab_ROM$ n_t1)
  weightedSD(data_lab_ROM$body_mass, data_lab_ROM$n_t1)
  
  
  ###########################################
  #Mean models - RR models
  ###########################################
  
  # Overall effect of temperature increase on dive duration
  
  model1.RR <- rma.mv(yi = yi, V = V, 
                      random = list(~1|study_ID, ~1|species_rotl, ~1|obs), 
                      R = list(species_rotl = PhyloA), test = "t",
                      data = data_lab_ROM)
  summary(model1.RR)
  
  I2(model1.RR, v = data_lab_ROM$vi, phylo = "species_rotl")

  
  # Model with moderators
  model2.RR <- rma.mv(yi = yi, V = V, 
                      mods = ~ mean_t + delta_t + log(body_mass) + respiration_mode, 
                      random = list(~1|study_ID, ~1|species_rotl, ~1|obs), 
                      R = list(species_rotl = PhyloA), test = "t",
                      data = data_lab_ROM)
  
 
  summary(model2.RR)  
  
  
  # Publication bias (Funnels):
  res <- residuals(model2.RR)
  funnel(res, vi = data_lab_ROM$vi, yaxis = "seinv")
  
  #Publication bias check- sampling variance (vi) included as a moderator
  model2.RR.pubbias <- rma.mv(yi = yi, V = V, 
                      mods = ~ mean_t + delta_t + log(body_mass) + respiration_mode + vi, 
                      random = list(~1|study_ID, ~1|species_rotl, ~1|obs), 
                      R = list(species_rotl = PhyloA), test = "t",
                      data = data_lab_ROM)
  summary(model2.RR.pubbias)
  
 
  # Do some model checks
  hist(residuals(model2.RR)) 
  plot(residuals(model2.RR))
  residuals(model2.RR)  

  
  
  ## Effect sizes for bimodal versus aerial breathers. Note here that you can fit the model with and without the intercept. Both results presented.
 
  data_lab_ROM %>% group_by(respiration_mode) %>% summarise(mean_t = mean(delta_t))
  data_lab_ROM %>% group_by(respiration_mode) %>% summarise(BM = mean(body_mass))
  data_lab_ROM %>% group_by(respiration_mode) %>% summarise(BM = sd(body_mass))
  
  model3.RR <- rma.mv(yi = yi, V = V, 
                      mods = ~ mean_t + delta_t + log(body_mass) + respiration_mode, 
                      random = list(~1|study_ID, ~1|species_rotl, ~1|obs), 
                      R = list(species_rotl = PhyloA), test = "t", 
                      data = data_lab_ROM)
  summary(model3.RR)
  
  # Do some model checks
  hist(residuals(model3.RR)) 
  plot(residuals(model3.RR))
  residuals(model3.RR)
  
  
  ## Effect sizes for different magnitudes of temperature increase (+3, +5-7, +8-9, +>10C). Controlling for the average temperature of the treatments
  
  model4.RR <- rma.mv(yi = yi, V = V, 
                      mods = ~ mean_t + t_magnitude-1, 
                      random = list(~1|study_ID, ~1|species_rotl, ~1|obs), 
                      R = list(species_rotl = PhyloA), test = "t",
                      data = data_lab_ROM)
  summary(model4.RR)
  
  # Model check
  hist(residuals(model4.RR)) 
  

  
  ###########################################
  #Variance models – CVR models
  ###########################################
  
  #Overall effect of temperature increase on dive duration variability
  model5.CVR <- rma.mv(yi = yi, V = V2, 
                       random = list(~1|study_ID, ~1|species_rotl, ~1|obs), 
                       R = list(species_rotl = PhyloA), test = "t",
                       data = data_lab_CVR)
  summary(model5.CVR)
  
  I2(model5.CVR, v = data_lab_CVR$vi, phylo = "species_rotl")  
  
  #Model with moderators
  model6.CVR <- rma.mv(yi = yi, V = V2, 
                       mods = ~ mean_t + delta_t + log(body_mass) + respiration_mode, 
                       random = list(~1|study_ID, ~1|species_rotl, ~1|obs), 
                       R = list(species_rotl = PhyloA), test = "t", 
                       data = data_lab_CVR)
  summary(model6.CVR)
  
  # Do some model checks
  hist(residuals(model6.CVR))
  plot(residuals(model6.CVR))
  residuals(model6.CVR)
  
  # Publication bias (Funnels):
  res <- residuals(model6.CVR)
  #funnel(model2.RR, level = c(0.90, 0.95, 0.99), yaxis = "seinv")
  funnel(res, vi = data_lab_CVR$vi, yaxis = "seinv")
  
  #Publication bias check- sampling variance (vi) included as a moderator
  model6.CVR.pubbias <- rma.mv(yi = yi, V = V2, 
                       mods = ~ mean_t + delta_t + log(body_mass) + respiration_mode + vi, 
                       random = list(~1|study_ID, ~1|species_rotl, ~1|obs), 
                       R = list(species_rotl = PhyloA), test = "t", 
                       data = data_lab_CVR)
  summary(model6.CVR.pubbias)
  
  
 
  
  ##Effect of respiration mode on CVR
  model7.CVR <- rma.mv(yi = yi, V = V2, 
                       mods = ~ mean_t + delta_t + log(body_mass) + respiration_mode-1, 
                       random = list(~1|study_ID, ~1|species_rotl, ~1|obs), 
                       R = list(species_rotl = PhyloA), test = "t", 
                       data = data_lab_CVR)
  
  model7.CVR <- rma.mv(yi = yi, V = V2, 
                       mods = ~ mean_t + delta_t + log(body_mass) + respiration_mode, 
                       random = list(~1|study_ID, ~1|species_rotl, ~1|obs), 
                       R = list(species_rotl = PhyloA), test = "t", 
                       data = data_lab_CVR)
  summary(model7.CVR)
  
  # Do some model checks
  hist(residuals(model7.CVR))

  
  ##Effect sizes for different magnitudes of temperature increase (+3, +5-7, +8-9, +>10C)
  model8.CVR <- rma.mv(yi = yi, V = V2, 
                       mods = ~ mean_t + t_magnitude-1, 
                       random = list(~1|study_ID, ~1|species_rotl, ~1|obs), 
                       R = list(species_rotl = PhyloA), test = "t", 
                       data = data_full_CVR)
  summary(model8.CVR)
  
  # Do some model checks
  hist(residuals(model8.CVR))
    