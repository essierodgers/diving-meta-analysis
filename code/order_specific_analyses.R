# Packages
pacman::p_load(metafor, MCMCglmm, tidyverse, rotl, phytools, corrplot, ape)

## Read some data
data <- read.csv("./data/lab_dive_durations_armbased.csv")
data <- read.csv("~/diving-meta-analysis/data/order_data/serpentes.csv")

#importing phylogeny from timetree
tree <- read.tree("~/diving-meta-analysis/data/order_data/serpentes_phylogeny.NWK")
plot(tree)
mat <- MCMCglmm::inverseA(tree, nodes="TIPS")$Ainv

# Separating genus and species
genus <- as.data.frame(do.call("rbind", str_split(str_trim(data$species, side = "both"), " ")))
names(genus) <- c("genus", "species_new")
data <- cbind(data, genus)	


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



# Bayesian with Tcentering                
data2 <- data[complete.cases(data[,c("acclimation_temp", "body_mass_g")]),]
prior_slope <- list(R = list(V = 1, nu = 0.001),
                    G = list(G1 = list(V=1, nu = 0.02),
                             G2 = list(V = diag(2), nu = 2)))

prior_int <- list(R = list(V = 1, nu = 0.001),
                  G = list(G1 = list(V=1, nu = 0.02),
                           G2 = list(V = 1, nu = 0.02)))







########## Phylogeny
### Access taxon relationships from Open Tree of Life
# Match species in dataset
data$species_rotl <- paste0(data$genus, "_", data$species_new)
data2$species_rotl <- paste0(data2$genus, "_", data2$species_new)

tree <- tnrs_match_names(names = unique(data2$species_rotl), context_name = "Animals")
write.csv(unique(data2$species_rotl), file = "serpentes_species.csv")


#Create a tree based on itt id's found on the open tree of life
tl <- tol_induced_subtree(ott_ids=na.omit(tree$ott_id)) 
plot(tl)

#branch lengths with phytools
phylo_branch <-compute.brlen(tl, method = "Grafen", power= 0.5)

phylo_branch$tip.label <- gsub("_ott.*", "", phylo_branch$tip.label)
is.ultrametric(phylo_branch)
plot(phylo_branch)

#phylo variance-covariance matrix
mat <- MCMCglmm::inverseA(tree, "TIPS")$invA




data2$species_rotl <- ifelse(data2$species_rotl == "Natrix_maura", "Coluber_maurus", data2$species_rotl)



#Models with Tcentering- simple model
model1 <- rma.mv(log(mean) ~  scale(acclimation_temp) + scale(log(body_mass_g)) + T_w + T_b, V = mean_sv, random = list(~1 |study_ID, ~1|species, ~1|obs), data = data)
model2 <- rma.mv(log(mean) ~  scale(acclimation_temp) + scale(log(body_mass_g)) + T_w, V = mean_sv, random = list(~1 |study_ID, ~1|species, ~1|obs), data = data)







myTree <-ape::read.tree(text='(((Hydrophis_curtus:7.38740000,Hydrophis_elegans:7.38740000):44.00570600, (Natrix_maura:51.39310600):24.20689400, (Acrochordus_arafurae:75.60000000);')
