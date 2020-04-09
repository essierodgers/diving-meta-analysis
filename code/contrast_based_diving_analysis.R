
# Clean workspace
  rm(list = ls()) 

# Load libraries
  pacman::p_load(metafor, MCMCglmm, tidyverse, rotl, phytools, corrplot)
  library(devtools)
  install_github("daniel1noble/metaAidR"); library(metaAidR)
  source("./code/func.R")

#data processing
  rerun_data = FALSE
  if(rerun_data == TRUE){

    # Bring in the data and convert key variables to required classes. 
                  data <- read.csv("./data/lab_dive_durations_contrast.csv", stringsAsFactors = FALSE)
      data$body_mass_g <- as.numeric(data$body_mass_g)
      data$t_magnitude <- ordered(data$t_magnitude, levels = c("3", "5-7", "8-9", "10plus"))
         data$study_ID <- as.factor(data$study_ID)
    
    # Separate verts from inverts
      data_verts <- data %>% filter(study_ID != 13)
    
    # Separate genus and species
              genus <- as.data.frame(do.call("rbind", str_split(str_trim(data_verts$species, side = "both"), " ")))
       names(genus) <- c("genus", "species_new", "sub_spp")
         data_verts <- cbind(data_verts, genus)
    
    # Fix up the species names so they match with phylogeny
      data_verts$species_rotl <- paste0(data_verts$genus, "_", data_verts$species_new)
      data_verts <- data_verts[,-31]

    # Write the file, so it can be loaded more easily 
        write.csv(data_verts, "./data/data_verts.csv", row.names=FALSE)
  }else {
                  data_verts <- read.csv("./data/data_verts.csv")
      data_verts$t_magnitude <- ordered(data_verts$t_magnitude, 
                                  levels = c("3", "5-7", "8-9", "10plus"))
  }

 # Import TimeTree phylogeny
        tree <- read.tree("./data/order_data/vert_phylogeny.NWK")
        plot(tree)

    # Create a phylogenetic correlation matrix
        PhyloA <- vcv(tree, corr = TRUE)
    
    # Fix what is different in data
    data_verts$species_rotl <- ifelse(data_verts$species_rotl == "Chrysemys_dorbignyi", "Trachemys_dorbigni", data_verts$species_rotl)
    data_verts$species_rotl <- ifelse(data_verts$species_rotl == "Triturus_alpestris", "Ichthyosaura_alpestris", data_verts$species_rotl)
    

# Have a look at the mean-variance relationship
  plot_func(data_verts, "sd_control", "mean_control")

# Calculating effect sizes, measure="ROM" means lnRR es, append=true means will add es to data set. yi= es estimate, vi sampling variance.
  data_verts_ROM <- escalc(m1i = mean_treatment, 
                           m2i = mean_control, 
                           n1i = n_treatment, 
                           n2i = n_control, 
                           sd1i = sd_treatment, 
                           sd2i = sd_treatment, 
                           append = TRUE, measure ="ROM", 
                           data = data_verts)

  data_verts_CVR <- escalc(m1i = mean_treatment, 
                           m2i = mean_control, 
                           n1i = n_treatment, 
                           n2i = n_control, 
                           sd1i = sd_treatment, 
                           sd2i = sd_treatment, 
                           append = TRUE, measure ="CVR", 
                           data = data_verts)

# Error calc-residual variance at the observation level as metafor does not add this by default
  data_verts_ROM$obs <- 1:dim(data_verts_ROM)[1]
  data_verts_CVR$obs <- 1:dim(data_verts_CVR)[1]

# Make shared control matrix
  #Shared control for ROM dataset
    data_verts_ROM$sc_cluster <- interaction(data_verts_ROM$study_ID, 
                                           data_verts_ROM$shared_control)

  # Create the SC V matrix
      V <- metaAidR::make_VCV_matrix(data_verts_ROM, "vi", "sc_cluster", 
                                  type = "vcv", rho = 0.5)
  
  # Look at the shared control (co)variance matrix
    corrplot(as.matrix(V), is.corr = FALSE)
    write.csv(V, file = "./output_matrices/sc_matrix_rom.csv")

  #Shared control for CVR dataset
    data_verts_CVR$sc_cluster <- interaction(data_verts_CVR$study_ID, 
                                              data_verts_CVR$shared_control)
  # Create the SC V matrix
    V2 <- make_VCV_matrix(data_verts_CVR, "vi", "sc_cluster", type = "vcv", rho = 0.5)
    
  # Look at the shared control (co)variance matrix
    corrplot(as.matrix(V), is.corr = FALSE)
    write.csv(V, file = "./output_matrices/sc_matrix_cvr.csv")

###########################################
#Mean models - RR models
###########################################

# Overall effect of temperature increase on dive duration
  model1_RR <- rma.mv(yi = yi, V = V, 
                   random = list(~1|study_ID, ~1|species_rotl, ~1|obs), 
                   R = list(species_rotl = PhyloA), test = "t",
                   data = data_verts_ROM)
  summary(model1_RR)

# Model with moderators
  model2_RR <- rma.mv(yi = yi, V = V, 
                   mods = ~ mean_t + delta_t + log(body_mass_g) + respiration_mode, 
                   random = list(~1|study_ID, ~1|species_rotl, ~1|obs), 
                   R = list(species_rotl = PhyloA), test = "t",
                   data = data_verts_ROM)
  summary(model2_RR)

  data_verts_ROM %>% group_by(respiration_mode) %>% summarise(mean_t = mean(mean_t))

  # Do some model checks
    hist(residuals(model2)) # outliers -3; maybe check this doesn't change anything.

# Effect sizes for bimodal versus aerial breathers
  model3_RR <- rma.mv(yi = yi, V = V, 
                   mods = ~mean_t + delta_t + respiration_mode-1, 
                   random = list(~1|study_ID, ~1|species_rotl, ~1|obs), 
                   R = list(species_rotl = PhyloA), test = "t", 
                   data = data_verts_ROM)
  summary(model3_RR)

# Effect sizes for different magnitudes of temperature increase (+3, +5-7, +8-9, +>10C). Controlling for the average temperature of the treatments
  model4_RR <- rma.mv(yi = yi, V = V, 
                   mods = ~ mean_t + t_magnitude-1, 
                   random = list(~1|study_ID, ~1|species_rotl, ~1|obs), 
                   R = list(species_rotl = PhyloA), test = "t",
                   data = data_verts_ROM)
  summary(model4_RR)
  # Model check
  hist(residuals(model4)) #  again, need to check that this -3 isn't messing with things

data_verts_ROM %>% group_by(t_magnitude) %>% summarise(mean_t = mean(scale(mean_t)))

# Order effect sizes
  model5_RR <- rma.mv(yi = yi, V = V, 
                   mods = ~mean_t + delta_t + order-1, 
                   random = list(~1|study_ID, ~1|species_rotl,  ~1|obs), 
                   R = list(species_rotl = PhyloA), test = "t",
                   data = data_verts_ROM)
  summary(model5_RR)

###########################################
#Variance models – CVR models
###########################################

#Overall effect of temperature increase on dive duration variability
  model1_CVR <- rma.mv(yi = yi, V = V2, 
                   random = list(~1|study_ID, ~1|species_rotl, ~1|obs), 
                   R = list(species_rotl = PhyloA), test = "t",
                   data = data_verts_CVR)
  summary(model1_CVR)

#Model with moderators
  model2_CVR <- rma.mv(yi = yi, V = V2, 
                   mods = ~ mean_t + delta_t + log(body_mass_g) + respiration_mode, 
                   random = list(~1|study_ID, ~1|species_rotl, ~1|obs), 
                   R = list(species_rotl = PhyloA), test = "t", 
                   data = data_verts_CVR)
  summary(model2_CVR)

#Effect sizes for bimodal versus aerial breathers
  model3_CVR <- rma.mv(yi = yi, V = V2, 
                   mods = ~respiration_mode-1, 
                   random = list(~1|study_ID, ~1|species_rotl, ~1|obs), 
                   R = list(species_rotl = PhyloA), test = "t", 
                   data = data_verts_CVR)
  summary(model3_CVR)

#Effect sizes for different magnitudes of temperature increase (+3, +5-7, +8-9, +>10C)
  model4_CVR <- rma.mv(yi = yi, V = V2, 
                   mods = ~mean_t + t_magnitude-1, 
                  random = list(~1|study_ID, ~1|species_rotl, ~1|obs), 
                  R = list(species_rotl = PhyloA), test = "t", 
                  data = data_verts_CVR)
  summary(model4_CVR)

#Order effect sizes
  model5_CVR <- rma.mv(yi = yi, V = V2, 
                   mods = ~order-1, 
                   random = list(~1|study_ID, ~1|species_rotl, ~1|obs), 
                   R = list(species_rotl = PhyloA), test = "t", 
                   data = data_verts_CVR)
  summary(model5_CVR)


### PLots

  devtools::install_github("itchyshin/orchard_plot", subdir = "orchaRd", force = TRUE, build_vignettes = TRUE)
  library(orchaRd)
  source("./code/revised_orchard.R")

  table_results <- mod_results(model5_RR, mod_cat = "order", mod_cont=c("scale(mean_t)", "scale(delta_t)"))
  print(table_results)

  orchard_plot(table_results, mod = "order", xlab = "Order", angle=45) + 
  annotate(geom = "text", label = paste0("n = ",as.character(study$n)), x= 1, y = c(1:4)+0.3)

  study <- data_verts_CVR %>% group_by(order) %>% summarise(n = length(unique(study_ID)))

  orchard_plot(table_results, mod = "order", xlab = "Order", angle=45) + 
  annotate(geom = "text", label = paste0("n = ",as.character(study$n)), x= 1, y = c(1:4)+0.3)

str(model5_RR)

extract = function(model){
table = data.frame(est = model$ beta, ci.lb = model$ ci.lb, ci.ub =model$ ci.ub)
return(table)
}
extract(model5_RR)
