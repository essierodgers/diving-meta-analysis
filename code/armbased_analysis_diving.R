# Packages
pacman::p_load(metafor, MCMCglmm, tidyverse, rotl, phytools, corrplot, ape)

source("./code/func.R")

#Data processing
rerun = FALSE

if(rerun == TRUE){
  
  # Bring in the data and convert key variables to required classes. 
  data <- read.csv("./data/lab_dive_durations_armbased.csv", stringsAsFactors = FALSE)
  data$body_mass_g <- as.numeric(data$body_mass_g)
  data$study_ID <- as.factor(data$study_ID)
  
  # Separate verts from inverts
  data_verts <- data %>% filter(study_ID != 13)
  
  # Separating genus and species
  genus <- as.data.frame(do.call("rbind", str_split(str_trim(data_verts$species, side = "both"), " ")))
  names(genus) <- c("genus", "species_new")
  data_verts <- cbind(data_verts, genus)
  data_verts$species_rotl <- paste0(data_verts$genus, "_", data_verts$species_new)
  data_verts <- data_verts[,-22]	
  
  
  # Tcentering 
  #Tw-centering of tempreature within each species across all studies on that species
  sp <- split(data_verts, data_verts$species)
  data_verts$T_w <- do.call("rbind", lapply(sp, function(x) scale(x$T, scale = FALSE)))
  
  # Add in observation level random effects
  data_verts$obs <- 1:dim(data_verts)[1]
  
   # Write the file, so it can be loaded more easily 
  write.csv(data_verts, "./data/data_verts_arm.csv", row.names=FALSE)

}else {
  
  data_verts <- read.csv("./data/data_verts_arm.csv")
  
}

 
#Importing phylogeny from TimeTree
    tree <- read.tree("./data/order_data/vert_phylogeny.NWK")
    plot(tree)
    A <- inverseA(tree, nodes = "TIPS")$Ainv
   
  # Fix up the species names so they match with phylogeny
      data_verts$species_rotl <- paste0(data_verts$genus, "_", data_verts$species_new)

  # Fix what is different in data
    data_verts$species_rotl <- ifelse(data_verts$species_rotl == "Chrysemys_dorbignyi", "Trachemys_dorbigni", data_verts$species_rotl)
    data_verts$species_rotl <- ifelse(data_verts$species_rotl == "Triturus_alpestris", "Ichthyosaura_alpestris", data_verts$species_rotl)
   
  
  # Check what is different-- different number of species (16 compared to 13)
      setdiff(unique(data_verts$species_rotl), sort(rownames(A)))
  

# Calculate the sampling variance for all mean estimates
  data_verts$mean_sv <- with(data_verts, m_sv(mean, sd, n))
  data_verts$sd_sv <- with(data_verts, sd_sv(n))

# Exploratory analysis 
  with(data_verts, plot(log(mean) ~ log(sd)))
  with(data_verts, hist(log(mean)))
  with(data_verts, hist(log(sd)))
  ggplot(data_verts, aes(y = log(mean), x=T_w, colour= species))+
    geom_line()

    data_verts %>% filter(mean >100)
    data_verts %>% filter(species_rotl == "Crocodylus_porosus")

# Have a look at some details about the temperatures
    Temps <- data2 %>% group_by (species_rotl) %>% summarise(Num_T = length(unique(T)), rangeL = range(T)[1], rangeU = range(T)[2], diff = range(T)[2] - range(T)[1])
      range(Temps$diff)

# Bayesian Priors              
data2 <- data_verts[complete.cases(data_verts[,c("T_w", "body_mass_g", "respiration_mode")]),]

prior_slope <- list(R = list(V = 1, nu = 0.001),
                    G = list(G1 = list(V=1, nu = 0.02),
                             G2 = list(V = diag(2), nu = 2)))

	###########################################
	#Mean models 
	###########################################	
	
	    model1.mean <- MCMCglmm(
                  log(mean) ~  T_w + log(body_mass_g) + respiration_mode*T_w, 
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
