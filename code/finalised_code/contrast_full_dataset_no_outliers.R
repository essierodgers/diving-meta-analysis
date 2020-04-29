# Clean workspace
rm(list = ls()) 


# Load libraries
pacman::p_load(metafor, MCMCglmm, tidyverse, rotl, phytools, corrplot, devtools)
install_github("daniel1noble/metaAidR"); library(metaAidR)
source("./code/func.R")
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
  data <- read.csv("./data/contrast_full_dataset.csv", stringsAsFactors = FALSE)
  data$body_mass <- as.numeric(data$body_mass)
  data$t_magnitude <- ordered(data$t_magnitude, levels = c("plus3", "plus5-7", "plus8-9", "plus10"))
  data$study_ID <- as.factor(data$study_ID)
  data$study_type <- as.factor(data$study_type)
  
  # Separate genus and species
  genus <- as.data.frame(do.call("rbind", str_split(str_trim(data$species, side = "both"), " ")))
  names(genus) <- c("genus", "species_new", "sub_spp")
  data <- cbind(data, genus)
  
  # Fix up the species names so they match with phylogeny
  data$species_rotl <- paste0(data$genus, "_", data$species_new)
  
  
  # Write the file, so it can be loaded more easily 
  write.csv(data, "./data/data_full.csv", row.names=FALSE)
}else {
  data_full <- read.csv("./data/data_full.csv", stringsAsFactors=FALSE)
  data_full$t_magnitude <- ordered(data_full$t_magnitude, 
                                    levels = c("plus3", "plus5-7", "plus8-9", "plus10"))
}




# Import TimeTree phylogeny
tree <- read.tree("./data/full_phylogeny.NWK")
plot(tree)



# Create a phylogenetic correlation matrix
PhyloA <- vcv(tree, corr = TRUE)



#Tree checks_full dataset
  setdiff(unique(data_full$species_rotl), sort(rownames(PhyloA)))
  # Fix what is different in data
  data_full$species_rotl <- ifelse(data_full$species_rotl == "Trachemys_dorbigni_Trachemys_dorbigni", "Trachemys_dorbigni", data_full$species_rotl)
  data_full$species_rotl <- ifelse(data_full$species_rotl == "Crocodylus_johnstoni", "Crocodylus_johnsoni", data_full$species_rotl)

  # Check same number of levels
  length(rownames(PhyloA))
  length(unique(data_full$species_rotl))
 



# Have a look at the mean-variance relationship
plot_func(data_full, "sd_t1", "mean_t1")

# Calculating effect sizes, measure="ROM" means lnRR es, append=true means will add es to data set. yi= es estimate, vi sampling variance.
data_full_ROM <- escalc(m1i = mean_t2, 
                         m2i = mean_t1, 
                         n1i = n_t2, 
                         n2i = n_t1, 
                         sd1i = sd_t2, 
                         sd2i = sd_t1, 
                         append = TRUE, measure ="ROM", 
                         data = data_full)

data_full_CVR <- escalc(m1i = mean_t2, 
                         m2i = mean_t1, 
                         n1i = n_t2, 
                         n2i = n_t1, 
                         sd1i = sd_t2, 
                         sd2i = sd_t1, 
                         append = TRUE, measure ="CVR", 
                         data = data_full)

# Error calc-residual variance at the observation level as metafor does not add this by default
data_full_ROM$obs <- 1:dim(data_full_ROM)[1]
data_full_CVR$obs <- 1:dim(data_full_CVR)[1]

# Make shared control matrix for ROM dataset
  data_full_ROM$sc_cluster <- interaction(data_full_ROM$study_ID, 
                                         data_full_ROM$shared_control)

  # Create the Shared control V matrix
  V <- metaAidR::make_VCV_matrix(data_full_ROM, "vi", "sc_cluster", 
                               type = "vcv", rho = 0.5)
  # Look at the shared control (co)variance matrix
  corrplot(as.matrix(V), is.corr = FALSE)
  write.csv(V, file = "./output_matrices/sc_matrix_rom.csv")


  # Make shared control matrix for CVR dataset
  data_full_CVR$sc_cluster <- interaction(data_full_CVR$study_ID, 
                                         data_full_CVR$shared_control)
  # Create the SC V matrix
  V2 <- make_VCV_matrix(data_full_CVR, "vi", "sc_cluster", type = "vcv", rho = 0.5)

  # Look at the shared control (co)variance matrix
  corrplot(as.matrix(V2), is.corr = FALSE)
  write.csv(V2, file = "./output_matrices/sc_matrix_cvr.csv")

  
  
  
  ###########################################
  #Mean models - RR models
  ###########################################
  
  # Overall effect of temperature increase on dive duration
  
  model1.RR <- rma.mv(yi = yi, V = V, 
                      random = list(~1|study_ID, ~1|species_rotl, ~1|obs), 
                      R = list(species_rotl = PhyloA), test = "t",
                      data = data_full_ROM)
  summary(model1.RR)
  
  I2(model1.RR, v = data_full_ROM$vi, phylo = "species_rotl")

  
  # Model with moderators
  model2.RR <- rma.mv(yi = yi, V = V, 
                      mods = ~ mean_t + delta_t + log(body_mass) + respiration_mode-1 + study_type, 
                      random = list(~1|study_ID, ~1|species_rotl, ~1|obs), 
                      R = list(species_rotl = PhyloA), test = "t",
                      data = data_full_ROM)
  
 
  summary(model2.RR)  
  
  
  # Publication bias (Funnels):
  res <- residuals(model2.RR)
  funnel(res, vi = data_full_ROM$vi, yaxis = "seinv")
  
  
  # Eggers regression: Significant intercept for meta-analytic residuals suggest publication bias if all sources of heterogeneity are accounted for.
  w <- 1 / data_full_ROM$vi # weight= Inverse sampling error; precision is the inverse of standard errors or the sqrt(w)
  o <- sqrt(w)*(res + data_full_ROM$vi)
  
  Egger <- lm(o ~ sqrt(w))
  summary(Egger)

  
  # Do some model checks
  hist(residuals(model2.RR)) 
  plot(residuals(model2.RR))
  residuals(model2.RR)  

  
  
  ## Effect sizes for bimodal versus aerial breathers. Note here that you can fit the model with and without the intercept. Both results presented.
 
  data_full_ROM %>% group_by(respiration_mode) %>% summarise(mean_t = weighted.mean(delta_t))
  weighted.mean(data_full_ROM$ delta_t, data_full_ROM$ n_t1)
  weighted.sd(data_full_ROM$ delta_t, data_full_ROM$ n_t1, na.rm = TRUE)
  wt.sd(data_full_ROM$delta_t, data_full_ROM$n_t1, na.rm = TRUE)
  
  model3.RR <- rma.mv(yi = yi, V = V, 
                      mods = ~ mean_t + delta_t + log(body_mass) + respiration_mode-1, 
                      random = list(~1|study_ID, ~1|species_rotl, ~1|obs), 
                      R = list(species_rotl = PhyloA), test = "t", 
                      data = data_full_ROM)
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
                      data = data_full_ROM)
  summary(model4.RR)
  
  # Model check
  hist(residuals(model4.RR)) 
  
  
  ## Effect sizes for lab versus field studies. 
  model4.RR <- rma.mv(yi = yi, V = V, 
                      mods = ~ mean_t + delta_t + body_mass + study_type-1, 
                      random = list(~1|study_ID, ~1|species_rotl, ~1|obs), 
                      R = list(species_rotl = PhyloA_outliers_removed), test = "t", 
                      data = data_full_ROM)
  model4.RR <- rma.mv(yi = yi, V = V_no_outlier, 
                      mods = ~ mean_t + delta_t + study_type, 
                      random = list(~1|study_ID, ~1|species_rotl, ~1|obs), 
                      R = list(species_rotl = PhyloA_outliers_removed), test = "t", 
                      data = data_full_ROM)
  summary(model4.RR)
  
  # Do some model checks
  hist(residuals(model4.RR)) 
  plot(residuals(model4.RR))
  residuals(model4.RR)

  
  ###########################################
  #Variance models – CVR models
  ###########################################
  
  #Overall effect of temperature increase on dive duration variability
  model5.CVR <- rma.mv(yi = yi, V = V2, 
                       random = list(~1|study_ID, ~1|species_rotl, ~1|obs), 
                       R = list(species_rotl = PhyloA), test = "t",
                       data = data_full_CVR)
  summary(model5.CVR)
  
  I2(model5.CVR, v = data_full_CVR$vi, phylo = "species_rotl")  
  
  #Model with moderators
  model6.CVR <- rma.mv(yi = yi, V = V2, 
                       mods = ~ mean_t + delta_t + log(body_mass) + respiration_mode, 
                       random = list(~1|study_ID, ~1|species_rotl, ~1|obs), 
                       R = list(species_rotl = PhyloA), test = "t", 
                       data = data_full_CVR)
  summary(model6.CVR)
  
  # Do some model checks
  hist(residuals(model6.CVR))
  plot(residuals(model6.CVR))
  residuals(model6.CVR)
  
  # Publication bias (Funnels):
  res <- residuals(model6.CVR)
  #funnel(model2.RR, level = c(0.90, 0.95, 0.99), yaxis = "seinv")
  funnel(res, vi = data_full_CVR$vi, yaxis = "seinv")
  
  
  # Eggers regression: Significant intercept for meta-analytic residuals suggest publication bias if all sources of heterogeneity are accounted for.
  w <- 1 / data_full_CVR$vi # weight= Inverse sampling error; precision is the inverse of standard errors or the sqrt(w)
  o <- sqrt(w)*(res + data_full_CVR$vi)
  
  Egger <- lm(o ~ sqrt(w))
  summary(Egger)  
  
  
  ##Effect of respiration mode on CVR
  model7.CVR <- rma.mv(yi = yi, V = V2, 
                       mods = ~ mean_t + delta_t + log(body_mass) + respiration_mode-1, 
                       random = list(~1|study_ID, ~1|species_rotl, ~1|obs), 
                       R = list(species_rotl = PhyloA), test = "t", 
                       data = data_full_CVR)
  
  model7.CVR <- rma.mv(yi = yi, V = V2, 
                       mods = ~ mean_t + delta_t + log(body_mass) + respiration_mode, 
                       random = list(~1|study_ID, ~1|species_rotl, ~1|obs), 
                       R = list(species_rotl = PhyloA), test = "t", 
                       data = data_full_CVR)
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
    