# Clean workspace
rm(list = ls()) 


# Load libraries
pacman::p_load(metafor, MCMCglmm, tidyverse, rotl, phytools, corrplot, devtools)
install_github("daniel1noble/metaAidR"); library(metaAidR)
source("./code/func.R")
library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library(stringr)
library(forcats)

  # Bring in the data and convert key variables to required classes. 
  data <- read.csv("./data/contrast_field_dataset.csv", stringsAsFactors = FALSE)
  data$body_mass <- as.numeric(data$body_mass)
  data$t_magnitude <- ordered(data$t_magnitude, levels = c("plus3", "plus5-7", "plus8-9", "plus10"))
  data$study_ID <- as.factor(data$study_ID)
  data$respiration_mode <- as.factor(data$respiration_mode)
  
  
  # Import TimeTree phylogeny
tree <- read.tree("./data/phylo_field.NWK")
plot(tree)

tree2 <- read.tree("./data/phylo_field_nooutlier.NWK")
plot(tree2)

# Create a phylogenetic correlation matrix
PhyloA <- vcv(tree, corr = TRUE)
PhyloA2 <- vcv(tree2, corr = TRUE)

#Tree checks_full dataset
setdiff(unique(data$species_rotl), sort(rownames(PhyloA)))
# Check same number of levels
length(rownames(PhyloA))
length(unique(data$species_rotl))



# Have a look at the mean-variance relationship
plot_func(data, "sd_t1", "mean_t1")

# Calculating effect sizes, measure="ROM" means lnRR es, append=true means will add es to data set. yi= es estimate, vi sampling variance.
data_ROM <- escalc(m1i = mean_t2, 
                       m2i = mean_t1, 
                       n1i = n_t2, 
                       n2i = n_t1, 
                       sd1i = sd_t2, 
                       sd2i = sd_t1, 
                       append = TRUE, measure ="ROM", 
                       data = data)

data_CVR <- escalc(m1i = mean_t2, 
                       m2i = mean_t1, 
                       n1i = n_t2, 
                       n2i = n_t1, 
                       sd1i = sd_t2, 
                       sd2i = sd_t1, 
                       append = TRUE, measure ="CVR", 
                       data = data)

# Error calc-residual variance at the observation level as metafor does not add this by defaul$obs <- 1:dim(data_field_ROM)[1]
data_ROM$obs <- 1:dim(data_ROM)[1]

data_CVR$obs <- 1:dim(data_CVR)[1]



# Make shared control matrix for ROM dataset
data_ROM$sc_cluster <- interaction(data_ROM$study_ID, 
                                       data_ROM$shared_control)

# Create the Shared control V matrix
V <- make_VCV_matrix(data_ROM, "vi", "sc_cluster", 
                               type = "vcv", rho = 0.5)


# Look at the shared control (co)variance matrix
corrplot(as.matrix(V), is.corr = FALSE)
write.csv(V, file = "./output_matrices/sc_matrix_rom_field.csv")


# Make shared control matrix for CVR dataset
data_CVR$sc_cluster <- interaction(data_CVR$study_ID, 
                                       data_CVR$shared_control)
# Create the SC V matrix
V2 <- make_VCV_matrix(data_CVR, "vi", "sc_cluster", type = "vcv", rho = 0.5)

# Look at the shared control (co)variance matrix
corrplot(as.matrix(V2), is.corr = FALSE)
write.csv(V2, file = "./output_matrices/sc_matrix_cvr_field.csv")


###########################################
#Mean models - RR models
###########################################

# Overall effect of temperature increase on dive duration

model1.RR <- rma.mv(yi = yi, V = V, 
                    random = list(~1|study_ID, ~1|species_rotl, ~1|obs), 
                    R = list(species_rotl = PhyloA), test = "t",
                    data = data_ROM)
summary(model1.RR)

I2(model1.RR, v = data_ROM$vi, phylo = "species_rotl")

# Do some model checks- two huge outliers identified
hist(residuals(model1.RR)) 
plot(residuals(model1.RR))
residuals(model1.RR)   

#remove outliers
data_ROM_outlier_removed <- data_ROM[-c(15,16),]

#update V matrix
V_no_outlier_ROM <- make_VCV_matrix(data_ROM_outlier_removed, "vi", "sc_cluster", 
                                type = "vcv", rho = 0.5)

#model without outliers
model1.RR.nooutlier <- rma.mv(yi = yi, V = V_no_outlier_ROM, 
                    random = list(~1|study_ID, ~1|species_rotl, ~1|obs), 
                    R = list(species_rotl = PhyloA2), test = "t",
                    data = data_ROM_outlier_removed)
summary(model1.RR.nooutlier)

# Do some model checks- looks much more sensible
hist(residuals(model1.RR.nooutlier)) 
plot(residuals(model1.RR.nooutlier))
   
#Cannot disentangle study and species effects
I2(model1.RR.nooutlier, v = data_ROM_outlier_removed$vi, phylo = "species_rotl")


# Publication bias (Funnels):
res <- residuals(model1.RR.nooutlier)
funnel(res, vi = data_ROM_outlier_removed$vi, yaxis = "seinv")

#Publication bias check- sampling variance (vi) included as a moderator
model2.RR.outlier_removed_pubbias <- rma.mv(yi = yi, V = V_no_outlier_ROM, 
                            mods = ~  mean_t + delta_t + vi, 
                            random = list(~1|study_ID, ~1|species_rotl, ~1|obs), 
                            R = list(species_rotl = PhyloA2), test = "t",
                            data = data_ROM_outlier_removed)
summary(model2.RR.outlier_removed_pubbias)



###########################################
#Variance models – CVR models
###########################################

#Overall effect of temperature increase on dive duration variability
model4.CVR <- rma.mv(yi = yi, V = V2, 
                     random = list(~1|study_ID, ~1|species_rotl, ~1|obs), 
                     R = list(species_rotl = PhyloA), test = "t",
                     data = data_CVR)
summary(model4.CVR)

#Cannot disentangle study effects from species effects
I2(model4.CVR, v = data_CVR$vi, phylo = "species_rotl")  


#Do some model checks
hist(residuals(model4.CVR)) 
plot(residuals(model4.CVR))

#Overall effect of temperature increase on dive duration variability- no outlier
#remove outliers
data_CVR_outlier_removed <- data_CVR[-c(15,16),]

#update V matrix
V_no_outlier_CVR <- make_VCV_matrix(data_CVR_outlier_removed, "vi", "sc_cluster", 
                                    type = "vcv", rho = 0.5)


model4.CVR.nooutlier <- rma.mv(yi = yi, V = V_no_outlier_CVR, 
                     random = list(~1|study_ID, ~1|species_rotl, ~1|obs), 
                     R = list(species_rotl = PhyloA2), test = "t",
                     data = data_CVR_outlier_removed)
summary(model4.CVR.nooutlier)


#Cannot disentangle study effects from species effects
I2(model4.CVR.nooutlier, v = data_CVR_outlier_removed$vi, phylo = "species_rotl")  


#Do some model checks
hist(residuals(model4.CVR)) 
plot(residuals(model4.CVR))


# Publication bias (Funnels):
res <- residuals(model4.CVR.nooutlier)
funnel(res, vi = data_CVR_outlier_removed$vi, yaxis = "seinv")



#Publication bias check- sampling variance (vi) included as a moderator
model4.CVR.outlier_removed_pubbias <- rma.mv(yi = yi, V = V_no_outlier_CVR, 
                                            mods = ~ vi, 
                                            random = list(~1|study_ID, ~1|species_rotl, ~1|obs), 
                                            R = list(species_rotl = PhyloA2), test = "t",
                                            data = data_CVR_outlier_removed)
summary(model4.CVR.outlier_removed_pubbias)

#Calculing weighted means and SD for delta_t- outliers removed
weighted.mean(data_ROM_outlier_removed$ delta_t, data_ROM_outlier_removed$ n_t1)
weightedSD(data_ROM_outlier_removed$delta_t, data_ROM_outlier_removed$n_t1)

#Calculing weighted means and SD for body mass- outliers removed
weighted.mean(data_ROM_outlier_removed$ body_mass, data_ROM_outlier_removed$ n_t1)
weightedSD(data_ROM_outlier_removed$body_mass, data_ROM_outlier_removed$n_t1)

