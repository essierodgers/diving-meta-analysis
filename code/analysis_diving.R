# Packages
pacman::p_load(metafor, MCMCglmm, tidyverse, rotl, phytools, corrplot, ape)


## Read some data
	data <- read.csv("./data/lab_dive_durations_armbased.csv")


#separate verts from inverts
	data_verts <- data %>% filter(study_ID != 13)

# Separating genus and species
	genus <- as.data.frame(do.call("rbind", str_split(str_trim(data_verts$species, side = "both"), " ")))
	names(genus) <- c("genus", "species_new")
	data_verts <- cbind(data_verts, genus)	
	data_verts <- data_verts[,-22]
	
	# Sampling variance for mean
	m_sv <- function(mean, sd, n){
		#page 148 Nakagawa et al 2015
		(sd^2) / (n*(mean^2))
	}

	sd_sv <- function(n){
		#page 148 Nakagawa et al 2015
		(1/2)*(n-1)
	}

# Now we can calculate the sampling variance for all mean estimates

	data_verts$mean_sv <- with(data_verts, m_sv(mean, sd, n))
	data_verts$sd_sv <- with(data_verts, sd_sv(n))

# Exploratory analysis 
	with(data_verts, plot(log(mean) ~ log(sd)))
	with(data_verts, hist(log(mean)))
	with(data_verts, hist(log(sd)))
	ggplot(data_verts, aes(y = log(mean), x=T, colour= species))+
	geom_point()

# Add in observation level random effects
	data_verts$obs <- 1:dim(data_verts)[1]

 #importing phylogeny from timetree
 tree <- read.tree("./data/order_data/vert_phylogeny.NWK")
 plot(tree)
 A <- inverseA(tree, nodes = "TIPS")$Ainv
 
#Tree checks function
  tree_checks <- function(data, tree, dataCol, type = c("checks", "prune")){
   type = match.arg(type)
   # How many unique species exist in data and tree
   Numbers <- matrix(nrow = 2, ncol = 1)
   Numbers[1,1] <- length(unique(data$spp_name_phylo)) 
   Numbers[2,1] <- length(tree$tip.label) 
   rownames(Numbers)<- c("Species in data:", "Species in tree:")
   # Missing species or species not spelt correct      
   species_list1= setdiff(sort(tree$tip.label), sort(unique(data[,dataCol])))
   species_list2= setdiff(sort(unique(data[,dataCol])), sort(tree$tip.label) )
   if(type == "checks"){
     return(list(SpeciesNumbers = data.frame(Numbers), 
                 Species_InTree_But_NotData=species_list1, 
                 Species_InData_But_NotTree=species_list2))
   }
   if(type == "prune"){
     if(length(species_list2) >=1) stop("Sorry, you can only prune a tree when you have no taxa existing in the data that are not in the tree")
     return(ape::drop.tip(tree, species_list1))
   }
 }
 
# Bayesian Priors              
  data2 <- data_verts[complete.cases(data_verts[,c("acclimation_temp", "body_mass_g")]),]
  prior_slope <- list(R = list(V = 1, nu = 0.001),
                      G = list(G1 = list(V=1, nu = 0.02),
                               G2 = list(V = diag(2), nu = 2)))
  
  prior_int <- list(R = list(V = 1, nu = 0.001),
                    G = list(G1 = list(V=1, nu = 0.02),
                             G2 = list(V = 1, nu = 0.02)))

  
  
#Adding species list
data_verts$species_list <- paste0(data_verts$genus, "_", data_verts$species_new)
data2$species_list <- paste0(data2$genus, "_", data2$species_new)
  

# Check same number of levels
length(rownames(A))
length(unique(data_verts$species_list))
length(unique(data2$species_list))

# Check what is different-- different number of species (16 compared to 13)
setdiff(unique(data2$species_list), sort(rownames(A)))

# Fix what is different species name
data2$species_list <- ifelse(data2$species_list == "Chrysemys_dorbignyi", "Trachemys_dorbigni", data2$species_list)

#check tree
tree_checks(data2, tree, dataCol = "species_list", type ="checks")


#Mean model 
model1 <- MCMCglmm(log(mean) ~  log(body_mass_g) + T, mev = data2$mean_sv, random = ~us(1):study_ID + us(1+T):species_list, ginverse = list(species_list = A), data = data2, prior = prior_slope, nitt = 50000, burnin = 10000, thin = 30)
summary(model1)
plot(model1)
  

#SD model
model2 <- MCMCglmm(log(sd) ~ log(mean) + log(body_mass_g) + T, mev = data2$sd_sv, random = ~us(1):study_ID + us(1+T):species_list, ginverse = list(species_list = A), data = data2, prior = prior_slope, nitt = 50000, burnin = 10000, thin = 30)
summary(model2)
plot(model2)
  
  
  
  

  




#OLD CODE- IGNORE FOR NOW

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
