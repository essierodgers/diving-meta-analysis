# Clean workspace
rm(list = ls()) 


# Load libraries
pacman::p_load(metafor, MCMCglmm, tidyverse, rotl, phytools, corrplot, devtools)
library(metaAidR)
source("./code/func.R")
library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library(stringr)
library(forcats)

  # Bring in the data and convert key variables to required classes. 
  data_field <- read.csv("./data/contrast_field_dataset.csv", stringsAsFactors = FALSE)
    data_field$body_mass <- as.numeric(data_field$body_mass)
  data_field$t_magnitude <- ordered(data_field$t_magnitude, 
                                    levels = c("plus3", "plus5-7", "plus8-9", "plus10"))
     data_field$study_ID <- as.factor(data_field$study_ID)
     data_field$respiration_mode <- as.factor(data_field$respiration_mode)
  
# Have a look at the mean-variance relationship
plot_func(data, "sd_t1", "mean_t1")

# Calculating effect sizes, measure="ROM" means lnRR es, append=true means will add es to data set. yi= es estimate, vi sampling variance.
data_field_ROM <- escalc(m1i = mean_t2, 
                       m2i = mean_t1, 
                       n1i = n_t2, 
                       n2i = n_t1, 
                       sd1i = sd_t2, 
                       sd2i = sd_t1, 
                       append = TRUE, measure ="ROM", 
                       data = data_field)

data_field_CVR <- escalc(m1i = mean_t2, 
                       m2i = mean_t1, 
                       n1i = n_t2, 
                       n2i = n_t1, 
                       sd1i = sd_t2, 
                       sd2i = sd_t1, 
                       append = TRUE, measure ="CVR", 
                       data = data_field)

# Error calc-residual variance at the observation level as metafor does not add this by defaul$obs <- 1:dim(data_field_ROM)[1]
data_field_ROM$obs <- 1:dim(data_field_ROM)[1]
data_field_CVR$obs <- 1:dim(data_field_CVR)[1]

# Make shared control matrix for ROM dataset
data_field_ROM$sc_cluster <- interaction(data_field_ROM$study_ID, 
                                       data_field_ROM$shared_control)

# Create the Shared control V matrix
V <- make_VCV_matrix(data_field_ROM, "vi", "sc_cluster", 
                               type = "vcv", rho = 0.5)


# Look at the shared control (co)variance matrix
corrplot(as.matrix(V), is.corr = FALSE)
write.csv(V, file = "./output_matrices/sc_matrix_rom_field.csv")


# Make shared control matrix for CVR dataset
data_field_CVR$sc_cluster <- interaction(data_field_CVR$study_ID, 
                                       data_field_CVR$shared_control)
# Create the SC V matrix
V2 <- make_VCV_matrix(data_field_CVR, "vi", "sc_cluster", type = "vcv", rho = 0.5)

# Look at the shared control (co)variance matrix
corrplot(as.matrix(V2), is.corr = FALSE)
write.csv(V2, file = "./output_matrices/sc_matrix_cvr_field.csv")


###########################################
#Mean models - RR models
###########################################

# Overall effect of temperature increase on dive duration

model1.field.RR <- rma.mv(yi = yi, V = V, 
                    random = list(~1|study_ID, ~1|obs), test = "t",
                    data = data_field_ROM)
summary(model1.field.RR)

# Do some model checks
hist(residuals(model1.field.RR)) 
plot(residuals(model1.field.RR))
residuals(model1.field.RR)   

# Publication bias (Funnels):
res <- residuals(model1.field.RR)
funnel(res, vi = model1.field.RR$vi, yaxis = "seinv")

#Publication bias check- sampling variance (vi) included as a moderator
model1.field.RR_pubias <- rma.mv(yi = yi, V = V, 
                            mods = ~  mean_t + delta_t + vi, 
                            random = list(~1|study_ID, ~1|obs), 
                             test = "t",
                            data = data_field_ROM)
summary(model1.field.RR_pubias)

###########################################
#Variance models – CVR models
###########################################

#Overall effect of temperature increase on dive duration variability
model4.field.CVR <- rma.mv(yi = yi, V = V2, 
                     random = list(~1|study_ID, ~1|obs), test = "t",
                     data = data_field_CVR)
summary(model4.field.CVR)

#Do some model checks
hist(residuals(model4.field.CVR)) 
plot(residuals(model4.field.CVR))

# Publication bias (Funnels):
res <- residuals(model4.field.CVR)
funnel(res, vi = model4.field.CVR$vi, yaxis = "seinv")

#Publication bias check- sampling variance (vi) included as a moderator
model4.field.CVR_pubbias <- rma.mv(yi = yi, V = V2, 
                                            mods = ~ vi, 
                                            random = list(~1|study_ID, ~1|obs), 
                                            test = "t",
                                            data = data_field_CVR)
summary(model4.field.CVR_pubbias)

#Calculating weighted means and SD for delta_t- outliers removed
weighted.mean(data_field_ROM$ delta_t, data_field_ROM$ n_t1)
weightedSD(data_field_ROM$delta_t, data_field_ROM$n_t1)

#Calculating weighted means and SD for body mass- outliers removed
weighted.mean(data_field_ROM$ body_mass, data_field_ROM$ n_t1)
weightedSD(data_field_ROM$body_mass, data_field_ROM$n_t1)

