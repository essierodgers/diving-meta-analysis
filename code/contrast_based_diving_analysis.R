
# Clean workspace
  rm(list = ls()) 

# Load libraries
  pacman::p_load(metafor, MCMCglmm, tidyverse, rotl, phytools, corrplot, devtools)

  install_github("daniel1noble/metaAidR"); library(metaAidR)
  source("./code/func.R")

#data processing
  rerun_data = FALSE
  if(rerun_data == TRUE){

    # Bring in the data and convert key variables to required classes. 
                  data <- read.csv("./data/lab_dive_durations_contrast.csv", stringsAsFactors = FALSE)
      data$body_mass_g <- as.numeric(data$body_mass_g)
      data$t_magnitude <- ordered(data$t_magnitude, levels = c("plus3", "plus5-7", "plus8-9", "plus10"))
         data$study_ID <- as.factor(data$study_ID)
    
    # Separate verts from inverts
      data_verts <- data %>% filter(study_ID != 13)
    
    # Separate genus and species
              genus <- as.data.frame(do.call("rbind", str_split(str_trim(data_verts$species, side = "both"), " ")))
       names(genus) <- c("genus", "species_new", "sub_spp")
         data_verts <- cbind(data_verts, genus)
    
    # Fix up the species names so they match with phylogeny
      data_verts$species_rotl <- paste0(data_verts$genus, "_", data_verts$species_new)
      

    # Write the file, so it can be loaded more easily 
        write.csv(data_verts, "./data/data_verts.csv", row.names=FALSE)
  }else {
                  data_verts <- read.csv("./data/data_verts.csv", stringsAsFactors=FALSE)
      data_verts$t_magnitude <- ordered(data_verts$t_magnitude, 
                                  levels = c("plus3", "plus5-7", "plus8-9", "plus10"))
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
  plot_func(data_verts, "sd_t1", "mean_t1")

# Calculating effect sizes, measure="ROM" means lnRR es, append=true means will add es to data set. yi= es estimate, vi sampling variance.
  data_verts_ROM <- escalc(m1i = mean_t2, 
                           m2i = mean_t1, 
                           n1i = n_t2, 
                           n2i = n_t1, 
                           sd1i = sd_t2, 
                           sd2i = sd_t1, 
                           append = TRUE, measure ="ROM", 
                           data = data_verts)

  data_verts_CVR <- escalc(m1i = mean_t2, 
                           m2i = mean_t1, 
                           n1i = n_t2, 
                           n2i = n_t1, 
                           sd1i = sd_t2, 
                           sd2i = sd_t1, 
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
      V <- make_VCV_matrix(data_verts_ROM, "vi", "sc_cluster", 
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

  model1.RR <- rma.mv(yi = yi, V = V, 
                   random = list(~1|study_ID, ~1|species_rotl, ~1|obs), 
                   R = list(species_rotl = PhyloA), test = "t",
                   data = data_verts_ROM)
  summary(model1.RR)

  I2(model1.RR, v = data_verts_ROM$vi, phylo = "species_rotl")

# Model with moderators
  model2.RR <- rma.mv(yi = yi, V = V, 
                   mods = ~ mean_t + delta_t + log(body_mass_g) + respiration_mode-1, 
                   random = list(~1|study_ID, ~1|species_rotl, ~1|obs), 
                   R = list(species_rotl = PhyloA), test = "t",
                   data = data_verts_ROM)
  
  data_verts_ROM %>% group_by(respiration_mode) %>% summarise(mean_t = mean(mean_t))

  summary(model2.RR)

  # Do some model checks
    hist(residuals(model2.RR)) # outliers -3; maybe check this doesn't change anything.
    plot(residuals(model2.RR))
    residuals(model2.RR)
  #remove outlier- need to updated V afterwards
    data_verts_ROM_outlier_removed <- data_verts_ROM[-48,]
 
   #update V matrix
    V_no_outlier <- make_VCV_matrix(data_verts_ROM_outlier_removed, "vi", "sc_cluster", 
                         type = "vcv", rho = 0.5)
    
    #rerun model without outlier
    model2.RR.nooutlier <- rma.mv(yi = yi, V = V_no_outlier, 
                        mods = ~ mean_t + delta_t + log(body_mass_g) + respiration_mode, 
                        random = list(~1|study_ID, ~1|species_rotl, ~1|obs), 
                        R = list(species_rotl = PhyloA), test = "t",
                        data = data_verts_ROM_outlier_removed)
    summary(model2.RR.nooutlier)
    #rerun model checks
    hist(residuals(model2.RR.nooutlier))
    
    
# Effect sizes for bimodal versus aerial breathers. Note here that you can fit the model with and without the intercept. Both results presented.
  model3.RR <- rma.mv(yi = yi, V = V, 
                   mods = ~ mean_t + delta_t + respiration_mode, 
                   random = list(~1|study_ID, ~1|species_rotl, ~1|obs), 
                   R = list(species_rotl = PhyloA), test = "t", 
                   data = data_verts_ROM)
  summary(model3.RR)

  # Do some model checks
  hist(residuals(model3.RR)) # outliers -3; maybe check this doesn't change anything.
  plot(residuals(model3.RR))
  residuals(model3.RR)
  
  #rerun without outlier
  model3.RR.nooutlier <- rma.mv(yi = yi, V = V_no_outlier, 
                      mods = ~ mean_t + delta_t + respiration_mode-1, 
                      random = list(~1|study_ID, ~1|species_rotl, ~1|obs), 
                      R = list(species_rotl = PhyloA), test = "t", 
                      data = data_verts_ROM_outlier_removed)
  summary(model3.RR.nooutlier)
  

# Effect sizes for different magnitudes of temperature increase (+3, +5-7, +8-9, +>10C). Controlling for the average temperature of the treatments

  model4.RR <- rma.mv(yi = yi, V = V, 
                   mods = ~ mean_t + t_magnitude-1, 
                   random = list(~1|study_ID, ~1|species_rotl, ~1|obs), 
                   R = list(species_rotl = PhyloA), test = "t",
                   data = data_verts_ROM)
  summary(model4.RR)
  
  # Model check
  hist(residuals(model4.RR)) #  again, need to check that this -3 isn't messing with things

  #rerun without outlier 
    model4.RR.nooutlier <- rma.mv(yi = yi, V = V_no_outlier, 
                      mods = ~ mean_t + t_magnitude-1, 
                      random = list(~1|study_ID, ~1|species_rotl, ~1|obs), 
                      R = list(species_rotl = PhyloA), test = "t",
                      data = data_verts_ROM_outlier_removed)
  summary(model4.RR.nooutlier)
  hist(residuals(model4.RR.nooutlier))
  
  
# Order effect sizes
  model5.RR <- rma.mv(yi = yi, V = V, 
                   mods = ~ mean_t + delta_t + order-1, 
                   random = list(~1|study_ID, ~1|species_rotl,  ~1|obs), 
                   R = list(species_rotl = PhyloA), test = "t",
                   data = data_verts_ROM)
  summary(model5.RR)
  
  # Model check
  hist(residuals(model5.RR)) #  again, need to check that this -3 isn't messing with things
  
  #rerun without outlier
  model5.RR.nooutlier <- rma.mv(yi = yi, V = V_no_outlier, 
                                mods = ~ mean_t + delta_t + order-1, 
                                random = list(~1|study_ID, ~1|species_rotl, ~1|obs), 
                                R = list(species_rotl = PhyloA), test = "t",
                                data = data_verts_ROM_outlier_removed)
  summary(model5.RR.nooutlier)
  hist(residuals(model5.RR.nooutlier))
  
###########################################
#Variance models – CVR models
###########################################

#Overall effect of temperature increase on dive duration variability
  model6.CVR <- rma.mv(yi = yi, V = V2, 
                   random = list(~1|study_ID, ~1|species_rotl, ~1|obs), 
                   R = list(species_rotl = PhyloA), test = "t",
                   data = data_verts_CVR)
  summary(model6.CVR)

I2(model6.CVR, v = data_verts_CVR$vi, phylo = "species_rotl")

#Model with moderators
  model7.CVR <- rma.mv(yi = yi, V = V2, 
                   mods = ~ mean_t + delta_t + log(body_mass_g) + respiration_mode-1, 
                   random = list(~1|study_ID, ~1|species_rotl, ~1|obs), 
                   R = list(species_rotl = PhyloA), test = "t", 
                   data = data_verts_CVR)
  summary(model7.CVR)

     # Do some model checks
      hist(residuals(model7.CVR)) # outlier around 3; maybe check this doesn't change anything.
      plot(residuals(model7.CVR))
      residuals(model7.CVR)
      
      #remove outlier- need to updated V afterwards
      data_verts_CVR_outlier_removed <- data_verts_CVR[-48,]
  
      #update V matrix
    V2_no_outlier <- make_VCV_matrix(data_verts_CVR_outlier_removed, 
                                    "vi", "sc_cluster", 
                                  type = "vcv", rho = 0.5)
  
      #rerun model without outlier
      model7.CVR.nooutlier <- rma.mv(yi = yi, V = V2_no_outlier, 
                           mods = ~ mean_t + delta_t + log(body_mass_g) + respiration_mode, 
                           random = list(~1|study_ID, ~1|species_rotl, ~1|obs), 
                           R = list(species_rotl = PhyloA), test = "t", 
                           data = data_verts_CVR_outlier_removed)
      summary(model7.CVR.nooutlier)
      hist(residuals(model7.CVR.nooutlier))
  
  model8.CVR <- rma.mv(yi = yi, V = V2, 
                   mods = ~ mean_t + delta_t + respiration_mode-1, 
                   random = list(~1|study_ID, ~1|species_rotl, ~1|obs), 
                   R = list(species_rotl = PhyloA), test = "t", 
                   data = data_verts_CVR)
  summary(model8.CVR)
  
  # Do some model checks
  hist(residuals(model8.CVR)) # outlier around 3; maybe check this doesn't change anything.
  
  #rerun model without outlier
  model8.CVR.nooutlier <- rma.mv(yi = yi, V = V2_no_outlier, 
                       mods = ~ mean_t + delta_t + respiration_mode-1, 
                       random = list(~1|study_ID, ~1|species_rotl, ~1|obs), 
                       R = list(species_rotl = PhyloA), test = "t", 
                       data = data_verts_CVR_outlier_removed)
  summary(model8.CVR.nooutlier)
  hist(residuals(model8.CVR.nooutlier))

#Effect sizes for different magnitudes of temperature increase (+3, +5-7, +8-9, +>10C)
  model9.CVR <- rma.mv(yi = yi, V = V2, 
                   mods = ~ mean_t + t_magnitude-1, 
                  random = list(~1|study_ID, ~1|species_rotl, ~1|obs), 
                  R = list(species_rotl = PhyloA), test = "t", 
                  data = data_verts_CVR)
  summary(model9.CVR)

  # Do some model checks
  hist(residuals(model9.CVR)) # outlier around 3; maybe check this doesn't change anything.
  
  
  #rerun model without outlier
    model9.CVR.nooutlier <- rma.mv(yi = yi, V = V2_no_outlier, 
                       mods = ~ mean_t + t_magnitude-1, 
                       random = list(~1|study_ID, ~1|species_rotl, ~1|obs), 
                       R = list(species_rotl = PhyloA), test = "t", 
                       data = data_verts_CVR_outlier_removed) 
  summary(model9.CVR.nooutlier) #only 10+ magnitude significant now
  hist(residuals(model9.CVR.nooutlier))
  
  
  #Order effect sizes
      model10.CVR <- rma.mv(yi = yi, V = V2, 
                   mods = ~ mean_t + delta_t + order-1, 
                   random = list(~1|study_ID, ~1|species_rotl, ~1|obs), 
                   R = list(species_rotl = PhyloA), test = "t", 
                   data = data_verts_CVR)
      summary(model10.CVR)

  
      # Do some model checks
      hist(residuals(model10.CVR)) # outlier around 3; maybe check this doesn't change anything.
      
      #rerun model without outlier
      model10.CVR.nooutlier <- rma.mv(yi = yi, V = V2_no_outlier, 
                            mods = ~ mean_t + delta_t + order-1, 
                            random = list(~1|study_ID, ~1|species_rotl, ~1|obs), 
                            R = list(species_rotl = PhyloA), test = "t", 
                            data = data_verts_CVR_outlier_removed)
      summary(model10.CVR.nooutlier)
      hist(residuals(model10.CVR.nooutlier))


###################################
## Figures
##################################

  devtools::install_github("itchyshin/orchard_plot", subdir = "orchaRd", force = TRUE, build_vignettes = TRUE)
  library(orchaRd)
  library(patchwork)
  source("./code/revised_orchard.R")


## Figure 1 mean and variance

  # Temperature moderator
  model1.RR_table_results <- mod_results(model1.RR, mod = "Int")
  print(model1.RR_table_results)
  model4.RR_table_results <- mod_results_new(model4.RR, mod_cat = "t_magnitude", mod_cont=c("mean_t"), type = "zero")
  print(model4.RR_table_results)

  model6.CVR_table_results <- mod_results(model6.CVR, mod = "Int")
  print(model6.CVR_table_results)
  model9.CVR_table_results <- mod_results_new(model9.CVR, mod_cat = "t_magnitude", mod_cont=c("mean_t"), type = "zero")
  print(model9.CVR_table_results)

  spp <- data_verts_ROM %>% group_by(t_magnitude) %>% summarise(n = length(unique(species_rotl)))

# Interesting issues here that none of us anticipated when making orchaRd and that is with respect to ordered factors. Things can get re-arranged in tables when order is not maintained. So, need to watch this. Fixed here, but colours are off. Just edit in Adobe

sppTotal <- length(unique(data_verts_ROM$species_rotl))
  p1_RR_mod1 <- orchard_plot(model1.RR_table_results, mod = "Int", xlab = "log Response Ratio (lnRR)", angle=45) + 
  annotate(geom = "text", label = paste0("italic(Sp) == ", sppTotal), x= 0.8, y = 1.295, size = 3.5, parse= TRUE) 
  p1_RR_mod4 <- orchard_plot(model4.RR_table_results, mod = "t_magnitude", xlab = "log Response Ratio (lnRR)", angle=45) + 
  annotate(geom = "text", label = paste0("italic(Sp) == ",as.character(spp$n)), x= 1, y = c(1:4)+0.29, size = 3.5, parse= TRUE) 

  p1_CVR_mod6 <- orchard_plot(model6.CVR_table_results, mod = "Int", xlab = "log Coefficient of Variation (lnCVR)", angle=45) + 
  annotate(geom = "text", label = paste0("italic(Sp) == ", sppTotal), x= 0.8, y = 1.295, size = 3.5, parse= TRUE) 
  p1_CVR_mod9 <- orchard_plot(model9.CVR_table_results, mod = "t_magnitude", xlab = "log Coefficient of Variation (lnCVR)", angle=45) + 
  annotate(geom = "text", label = paste0("italic(Sp) == ",as.character(spp$n)), x= 1, y = c(1:4)+0.29, size = 3.5, parse= TRUE) 

# Figure should be 75% there. Just do some modifications in Adobe.
  pdf(width = 12, height = 11, file = "./preliminary_figures/Fig1.pdf", useDingbats = FALSE)
    (p1_RR_mod1 / p1_RR_mod4) | (p1_CVR_mod6 / p1_CVR_mod9)
  dev.off()

################################
  ## Figure 2
################################
model2.RR
  model7.CVR

  model2.RR_table_results <- mod_results_new(model2.RR, mod_cat = "respiration_mode", mod_cont=c("mean_t", "delta_t", "log(body_mass_g)"), type = "mean")
  print(model2.RR_table_results)
  
  model7.CVR_table_results <- mod_results_new(model7.CVR, mod_cat = "respiration_mode", mod_cont=c("mean_t", "delta_t", "log(body_mass_g)"), type = "mean")
  print(model7.CVR_table_results)

  spp <- data_verts_ROM %>% group_by(respiration_mode) %>% summarise(n = length(unique(species_rotl)))

  p1_RR_mod2 <- orchard_plot(model2.RR_table_results, mod = "respiration_mode", xlab = "log Response Ratio (lnRR)", angle=45) + 
  annotate(geom = "text", label = paste0("italic(Sp) == ",as.character(spp$n)), x= 1, y = c(1:2)+0.29, size = 3.5, parse= TRUE) 

  p1_CVR_mod7 <- orchard_plot(model7.CVR_table_results, mod = "respiration_mode", xlab = "log Coefficient of Variation (lnCVR)", angle=45) + 
  annotate(geom = "text", label = paste0("italic(Sp) == ",as.character(spp$n)), x= 1, y = c(1:2)+0.29, size = 3.5, parse= TRUE)

pdf(width = 9.480176,  height = 4.546255, "./preliminary_figures/Fig2.pdf", useDingbats = FALSE)
  p1_RR_mod2 | p1_CVR_mod7
dev.off()