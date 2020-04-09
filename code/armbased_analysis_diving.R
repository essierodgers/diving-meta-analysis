# Packages
pacman::p_load(metafor, MCMCglmm, tidyverse, rotl, phytools, corrplot, ape)

source("./code/func.R")

#Data processing
rerun = TRUE
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
  
  
  #importing phylogeny from timetree
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
  
   # Write the file, so it can be loaded more easily 
  write.csv(data_verts, "./data/data_verts.csv", row.names=FALSE)
}else {
  data_verts <- read.csv("./data/data_verts.csv")
  
}


# Sampling variance for mean
m_sv <- function(mean, sd, n){
  #page 148 Nakagawa et al 2015
  (sd^2) / (n*(mean^2))
}

sd_sv <- function(n){
  #page 148 Nakagawa et al 2015
  (1/2)*(n-1)
}

# Calculate the sampling variance for all mean estimates
data_verts$mean_sv <- with(data_verts, m_sv(mean, sd, n))
data_verts$sd_sv <- with(data_verts, sd_sv(n))


# Exploratory analysis 
with(data_verts, plot(log(mean) ~ log(sd)))
with(data_verts, hist(log(mean)))
with(data_verts, hist(log(sd)))
ggplot(data_verts, aes(y = log(mean), x=T, colour= species))+
  geom_point()

# Bayesian Priors              
data2 <- data_verts[complete.cases(data_verts[,c("T_w", "body_mass_g", "respiration_mode")]),]
prior_slope <- list(R = list(V = 1, nu = 0.001),
                    G = list(G1 = list(V=1, nu = 0.02),
                             G2 = list(V = diag(2), nu = 2)))

prior_int <- list(R = list(V = 1, nu = 0.001),
                  G = list(G1 = list(V=1, nu = 0.02),
                           G2 = list(V = 1, nu = 0.02)))

	# Check same number of levels
	length(rownames(A))
	length(unique(data_verts$species_rotl))
	length(unique(data2$species_rotl))
	

	#check tree--not working anymore 
	tree_checks(data_verts, A, dataCol = "species_rotl", type ="checks")
	tree_checks(data2, A, dataCol = "species_rotl", type ="checks")
	
	
	

	###########################################
	#Mean models 
	###########################################	
	
	    model1.mean <- MCMCglmm(log(mean) ~  T_w + log(body_mass_g) + respiration_mode, 
	                            mev = data2$mean_sv, random = ~us(1):study_ID + us(1+T):species_rotl, 
	                            ginverse = list(species_rotl = A), 
	                            data = data2, prior = prior_slope, nitt = 130000, burnin = 30000, thin = 50)
	    summary(model1.mean)
	    plot(model1.mean)
	
    	#rerun with outlier
    	data_verts_outlier_removed <- data_verts[-48,]
	
    	# Bayesian Priors- without outlier             
    	data2_outlier_removed <- data_verts[complete.cases(data_verts_outlier_removed[,c("T_w", "body_mass_g", "respiration_mode")]),]
	    prior_slope <- list(R = list(V = 1, nu = 0.001),
	                    G = list(G1 = list(V=1, nu = 0.02),
	                             G2 = list(V = diag(2), nu = 2)))
	
	    prior_int <- list(R = list(V = 1, nu = 0.001),
	                  G = list(G1 = list(V=1, nu = 0.02),
	                           G2 = list(V = 1, nu = 0.02)))
	
    	#rerun model- doesn't change estimates very much
    	model1.mean.nooutlier <- MCMCglmm(log(mean) ~  T_w + log(body_mass_g) + respiration_mode,
    	                                  mev = data2_outlier_removed$mean_sv, random = ~us(1):study_ID + us(1+T):species_rotl,
    	                                  ginverse = list(species_rotl = A),
    	                                  data = data2_outlier_removed, prior = prior_slope, nitt = 130000, burnin = 30000, thin = 50)
	    summary(model1.mean.nooutlier)
	    plot(model1.mean.nooutlier)
	
	###########################################
	#Variance models – CVR models
	###########################################
	
	    #SD model
    	model2.sd <- MCMCglmm(log(sd) ~ log(mean) + T_w + log(body_mass_g) + respiration_mode,
    	                      mev = data2$sd_sv, random = ~us(1):study_ID + us(1+T):species_rotl, 
    	                      ginverse = list(species_rotl = A), 
    	                      data = data2, prior = prior_slope, nitt = 130000, burnin = 30000, thin = 50)
	    summary(model2.sd)
	    plot(model2.sd)

     #rerun model without outlier
    	model2.sd.nooutlier <- MCMCglmm(log(sd) ~ log(mean) + T_w + log(body_mass_g) + respiration_mode, 
    	                                mev = data2_outlier_removed$sd_sv, random = ~us(1):study_ID + us(1+T):species_rotl, 
    	                                ginverse = list(species_rotl = A), 
    	                                data = data2_outlier_removed, prior = prior_slope, nitt = 130000, burnin = 30000, thin = 50)
	 summary(model2.sd.nooutlier)
	 plot(model2.sd.nooutlier)




 


  
  





  
  
  
 

  




#OLD CODE- IGNORE FOR NOW
	#Adding species list
	data_verts$species_list <- paste0(data_verts$genus, "_", data_verts$species_new)
	data2$species_list <- paste0(data2$genus, "_", data2$species_new)
	
# Check what is different-- different number of species (16 compared to 13)
setdiff(unique(data2$species_rotl), sort(rownames(A)))

# Fix what is different species name
data2$species_rotl <- ifelse(data2$species_rotl == "Chrysemys_dorbignyi", "Trachemys_dorbigni", data2$species_rotl)
# Check what is different-- different number of species (16 compared to 13)
setdiff(unique(data2$species_rotl), sort(rownames(A)))

# Fix what is different species name
data2$species_rotl <- ifelse(data2$species_rotl == "Chrysemys_dorbignyi", "Trachemys_dorbigni", data2$species_rotl)



# Tcentering 
#Tw-centering of tempreature within each species across all studies on that species
sp <- split(data_verts, data_verts$species)
data_verts$T_w <- do.call("rbind", lapply(sp, function(x) scale(x$T, scale = FALSE)))


#Models with Tcentering- simple model
model1 <- rma.mv(log(mean) ~  scale(acclimation_temp) + scale(log(body_mass_g)) + T_w + T_b, V = mean_sv, random = list(~1 |study_ID, ~1|species, ~1|obs), data = data_verts)

# Bayesian with Tcentering                
data2 <- data[complete.cases(data[,c("acclimation_temp", "body_mass_g")]),]
prior_slope <- list(R = list(V = 1, nu = 0.001),
                    G = list(G1 = list(V=1, nu = 0.02),
                             G2 = list(V = diag(2), nu = 2)))

prior_int <- list(R = list(V = 1, nu = 0.001),
                  G = list(G1 = list(V=1, nu = 0.02),
                           G2 = list(V = 1, nu = 0.02)))

model3 <- MCMCglmm(log(mean) ~ scale(acclimation_temp) + scale(log(body_mass_g)) + T_w + T_b, mev = data2$mean_sv, random = ~us(1):study_ID + us(1):species, data = data2, prior = prior_int, nitt = 50000, burnin = 10000, thin = 30)
summary(model3)
plot(model3)





model4 <- MCMCglmm(log(mean) ~ scale(acclimation_temp) + scale(log(body_mass_g)) + T_w + T_b, mev = data2$mean_sv, random = ~us(1):study_ID + us(1+T):species_rotl, ginverse = list(species_rotl = A), data = data2, prior = prior_slope, nitt = 50000, burnin = 10000, thin = 30)
summary(model4)

model5 <- MCMCglmm(log(mean) ~ scale(acclimation_temp) + scale(log(body_mass_g)) + T_w, mev = data2$mean_sv, random = ~us(1):study_ID + us(1+T):species_rotl, ginverse = list(species_rotl = A), data = data2, prior = prior_slope, nitt = 50000, burnin = 10000, thin = 30)
summary(model5)
saveRDS(model5, "model5.rds")

model6 <- MCMCglmm(log(mean) ~ scale(log(body_mass_g)) + T, mev = data2$mean_sv, random = ~us(1):study_ID + us(1+T):species_rotl, ginverse = list(species_rotl = A), data = data2, prior = prior_slope, nitt = 50000, burnin = 10000, thin = 30)
summary(model5)

#simple model-verts only                 
model1 <- rma.mv(log(mean) ~  scale(acclimation_temp) + scale(log(body_mass_g)) + T_w + T_b, V = mean_sv, random = list(~1 |study_ID, ~1|species, ~1|obs), data = data_verts)                

#Phlogeny- verts only
data_verts$species_rotl <- paste0(data_verts$genus, "_", data_verts$species_new)
tree <- tnrs_match_names(names = unique(data_verts$species_rotl), context_name = "Animals")
write.csv(unique(data_verts$species_rotl), file = "species.csv")
#Create a tree based on itt id's found on the open tree of life
tl_verts <- tol_induced_subtree(ott_ids=na.omit(tree$ott_id)) 
plot(tl_verts)
#branch lengths with phytools
phylo_branch <-compute.brlen(tl_verts, method = "Grafen", power= 0.8)
phylo_branch$tip.label <- gsub("_ott.*", "", phylo_branch$tip.label)
is.ultrametric(phylo_branch)
plot(phylo_branch)
#phylo variance-covariance matrix
A <- inverseA(phylo_branch, nodes = "TIPS")$Ainv
data_verts$species_rotl <- ifelse(data_verts$species_rotl == "Chrysemys_dorbignyi", "Trachemys_dorbigni", data_verts$species_rotl)

# Bayesian with Tcentering -verts only

data2 <- data_verts[complete.cases(data_verts[,c("acclimation_temp", "body_mass_g")]),]
prior_slope <- list(R = list(V = 1, nu = 0.001),
                    G = list(G1 = list(V=1, nu = 0.02),
                             G2 = list(V = diag(2), nu = 2)))

prior_int <- list(R = list(V = 1, nu = 0.001),
                  G = list(G1 = list(V=1, nu = 0.02),
                           G2 = list(V = 1, nu = 0.02)))


model3 <- MCMCglmm(log(mean) ~ scale(acclimation_temp) + scale(log(body_mass_g)) + T_w + T_b, mev = data_verts$mean_sv, random = ~us(1):study_ID + us(1):species, data = data_verts, prior = prior_int, nitt = 50000, burnin = 10000, thin = 30)
summary(model3)
plot(model3)

model4 <- MCMCglmm(log(mean) ~ scale(acclimation_temp) + scale(log(body_mass_g)) + T_w  + T_b + I(T_b^2), mev = data_verts$mean_sv, random = ~us(1):study_ID + us(1+T):species_rotl, ginverse = list(species_rotl = A), data = data_verts, prior = prior_slope, nitt = 50000, burnin = 10000, thin = 30)
summary(model4)




model5 <- MCMCglmm(log(mean) ~ scale(acclimation_temp) + scale(log(body_mass_g)) + T_w + T_b, mev = data2$mean_sv, random = ~us(1):study_ID + us(1+T):species_rotl, ginverse = list(species_rotl = A), data = data2, prior = prior_slope, nitt = 50000, burnin = 10000, thin = 30)
summary(model5)

model7 <- MCMCglmm(log(sd) ~ log(mean) + log(body_mass_g) + T, mev = arms_data2$sd_sv, random = ~us(1):study_ID + us(1+T):species_rotl, ginverse = list(species_rotl = A), data = arms_data2, prior = prior_slope, nitt = 50000, burnin = 10000, thin = 30)
summary(model7)





#WIHTOUT Tcentering
# Now lets try a simple model
model1 <- rma.mv(log(mean) ~ scale(acclimation_temp) + scale(body_mass_g) + T + I(T^2), V = mean_sv, random = list(~1 |study_ID, ~1|species, ~1|obs), data = data)
model2 <- rma.mv(log(mean) ~ scale(body_mass_g) + T + I(T^2), V = mean_sv, random = list(~1 |study_ID, ~1|species, ~1|obs), data = data)

# Bayesian
data2 <- data_verts[complete.cases(data_verts[,c("acclimation_temp", "body_mass_g")]),]
prior_slope <- list(R = list(V = 1, nu = 0.001),
                    G = list(G1 = list(V=1, nu = 0.02),
                             G2 = list(V = diag(2), nu = 2)))

prior_int <- list(R = list(V = 1, nu = 0.001),
                  G = list(G1 = list(V=1, nu = 0.02),
                           G2 = list(V = 1, nu = 0.02)))




model3 <- MCMCglmm(log(mean) ~ scale(acclimation_temp) + scale(body_mass_g) + T_w + T_b, mev = data2$mean_sv, random = ~us(1):study_ID + us(1):species, data = data2, prior = prior_int, nitt = 50000, burnin = 10000, thin = 30)
summary(model3)
plot(model3)

# This is the one to go with. ginverse = phylogeny
model3 <- MCMCglmm(log(mean) ~ scale(acclimation_temp) + scale(body_mass_g) + T, mev = data2$mean_sv, random = ~us(1):study_ID + us(1+T):species, data = data2, prior = prior_slope, nitt = 50000, burnin = 10000, thin = 30)
summary(model3)



model5 <- MCMCglmm(log(mean) ~ scale(acclimation_temp) + scale(body_mass_g) + T*respiration_mode, mev = data2$mean_sv, random = ~us(1):study_ID + us(1+T):species, data = data2, prior = prior_slope, nitt = 50000, burnin = 10000, thin = 30)
summary(model5)

model4 <- MCMCglmm(log(mean) ~ scale(acclimation_temp) + scale(body_mass_g) + T, mev = data2$mean_sv, random = ~us(1):study_ID + us(1):species, data = data2, prior = prior_int, nitt = 50000, burnin = 10000, thin = 30)
summary(model4)





# Check what is different-- different number of species (16 compared to 13)
setdiff(unique(data2$species_rotl), sort(rownames(PhyloA)))

########## Phylogeny
### Access taxon relationships from Open Tree of Life
# Match species in dataset
data$species_rotl <- paste0(data$genus, "_", data$species_new)
data2$species_rotl <- paste0(data2$genus, "_", data2$species_new)

tree <- tnrs_match_names(names = unique(data2$species_rotl), context_name = "Animals")
write.csv(unique(data2$species_rotl), file = "species.csv")


#Create a tree based on itt id's found on the open tree of life
tl <- tol_induced_subtree(ott_ids=na.omit(tree$ott_id)) 
plot(tl)

#branch lengths with phytools
phylo_branch <-compute.brlen(tl, method = "Grafen", power= 0.8)

phylo_branch$tip.label <- gsub("_ott.*", "", phylo_branch$tip.label)
is.ultrametric(phylo_branch)
plot(phylo_branch)




#phylo variance-covariance matrix
A <- inverseA(phylo_branch, nodes = "TIPS")$Ainv

data2$species_rotl <- ifelse(data2$species_rotl == "Chrysemys_dorbignyi", "Trachemys_dorbigni", data2$species_rotl)
data2$species_rotl <- ifelse(data2$species_rotl == "Ilybius_erichsonii", "Ilybius_erichsoni", data2$species_rotl)
data2$species_rotl <- ifelse(data2$species_rotl == "Ilybius_pederzani", "Ilybius_pederzanii", data2$species_rotl)

sort(unique(data2$species_rotl)) == sort(rownames(A))

sort(rownames(A))[-which(sort(unique(data2$species_rotl)) == sort(rownames(A)))]
sort(unique(data2$species_rotl))[-which(sort(unique(data2$species_rotl)) == sort(rownames(A)))]
