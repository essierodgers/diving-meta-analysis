# Packages
pacman::p_load(metafor, MCMCglmm, tidyverse, rotl, phytools, corrplot, ape)

# Remove stuff
rm(list=ls())

#separate verts from inverts
data_verts <- data2 %>% filter(study_ID != 13)

## Read some data
	data <- read.csv("./data/lab_dive_durations_armbased.csv")
	data <- read.csv("~/diving-meta-analysis/data/lab_dive_durations_armbased.csv")


# Separating genus and species
	genus <- as.data.frame(do.call("rbind", str_split(str_trim(data$species, side = "both"), " ")))
	names(genus) <- c("genus", "species_new")
	data <- cbind(data, genus)	
	data <- data[,-22]
	
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

	data$mean_sv <- with(data, m_sv(mean, sd, n))
	  data$sd_sv <- with(data, sd_sv(n))

# Exploratory analysis 
	with(data, plot(log(mean) ~ log(sd)))
	with(data, hist(log(mean)))
	with(data, hist(log(sd)))
	ggplot(data, aes(y = log(mean), x=T, colour= species))+
	geom_point()

# Add in observation level random effects
	data$obs <- 1:dim(data)[1]

# Tcentering 
#Tw-centering of tempreature within each species across all studies on that species
	sp <- split(data, data$species)
	data$T_w <- do.call("rbind", lapply(sp, function(x) scale(x$T, scale = FALSE)))


#Models with Tcentering- simple model
model1 <- rma.mv(log(mean) ~  scale(acclimation_temp) + scale(log(body_mass_g)) + T_w + T_b, V = mean_sv, random = list(~1 |study_ID, ~1|species, ~1|obs), data = data)

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

#simple model-verts only                 
model1 <- rma.mv(log(mean) ~  scale(acclimation_temp) + scale(log(body_mass_g)) + T_w + T_b, V = mean_sv, random = list(~1 |study_ID, ~1|species, ~1|obs), data = data_verts)                

#Phlogeny- verts only
data_verts$species_rotl <- paste0(data$genus, "_", data$species_new)
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

data2 <- data[complete.cases(data[,c("acclimation_temp", "body_mass_g")]),]
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




#WIHTOUT Tcentering
# Now lets try a simple model
model1 <- rma.mv(log(mean) ~ scale(acclimation_temp) + scale(body_mass_g) + T + I(T^2), V = mean_sv, random = list(~1 |study_ID, ~1|species, ~1|obs), data = data)
model2 <- rma.mv(log(mean) ~ scale(body_mass_g) + T + I(T^2), V = mean_sv, random = list(~1 |study_ID, ~1|species, ~1|obs), data = data)

# Bayesian
data2 <- data[complete.cases(data[,c("acclimation_temp", "body_mass_g")]),]
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



model5 <- MCMCglmm(log(mean) ~ scale(acclimation_temp) + scale(log(body_mass_g)) + T_w + T_b, mev = data2$mean_sv, random = ~us(1):study_ID + us(1+T):species_rotl, ginverse = list(species_rotl = A), data = data2, prior = prior_slope, nitt = 50000, burnin = 10000, thin = 30)
summary(model5)
