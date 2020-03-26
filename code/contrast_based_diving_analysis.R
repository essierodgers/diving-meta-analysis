# Load libraries
pacman::p_load(metafor, MCMCglmm, tidyverse, rotl, phytools, corrplot)
library(devtools)
install_github("daniel1noble/metaAidR"); library(metaAidR)

# Read in the data
data <- read.csv("./data/lab_dive_durations_contrast_highlow.csv")
str(data)
data$body_mass_g <- as.numeric(data$body_mass_g)

#separate verts from inverts
data_verts <- data %>% filter(study_ID != 13)


#how many species per study
dat <- data %>%
  group_by(study_ID) %>%
  summarise(n = length(unique(species))) 

#importing phylogeny from timetree
tree <- read.tree("./data/order_data/vert_phylogeny.NWK")
plot(tree)
PhyloA <- vcv(tree, corr = TRUE)

# Separating genus and species
genus <- as.data.frame(do.call("rbind", str_split(str_trim(data_verts$species, side = "both"), " ")))
names(genus) <- c("genus", "species_new")
data_verts <- cbind(data_verts, genus)

data_verts$species_rotl <- paste0(data_verts$genus, "_", data_verts$species_new)

# Check same number of levels
length(rownames(PhyloA))
length(unique(data_verts$species_rotl))

# Check what is different
setdiff(unique(data_verts$species_rotl), sort(rownames(PhyloA)))

# Fix what is different in data
data_verts$species_rotl <- ifelse(data_verts$species_rotl == "Chrysemys_dorbignyi", "Trachemys_dorbigni", data_verts$species_rotl)
data_verts$species_rotl <- ifelse(data_verts$species_rotl == "Triturus_alpestris", "Ichthyosaura_alpestris", data_verts$species_rotl)


## Plots
plot_func <- function(data, x, y){
  plot(log(data[,y]) ~log(data[,x]), pch = 16, ylab = "log(y)", xlab = "log(x)") 
}

plot_func(data, "sd_control", "mean_control")


#calculating effect sizes, measure="ROM" means lnRR es, append=true means will add es to data set. yi= es estimate, vi sampling variance.
data_verts_ROM <- escalc(m1i = mean_control, m2i = mean_treatment, n1i = n_control, n2i = n_treatment, sd1i = sd_control, sd2i = sd_treatment, append = TRUE, measure ="ROM", data = data_verts)
data_verts_CVR <- escalc(m1i = mean_control, m2i = mean_treatment, n1i = n_control, n2i = n_treatment, sd1i = sd_control, sd2i = sd_treatment, append = TRUE, measure ="CVR", data = data_verts)





# Check directionality – control on numerator
(log(data$mean_control) - log(data$mean_treatment) ) == (log(data$mean_control/data$mean_treatment))

#error calc-residual variance at the observation level
data_verts_ROM$obs <- 1:dim(data_verts_ROM)[1]

# Meta-analytic multivariate model - wihtout intercept
model1 <- rma.mv(yi ~ -1 + scale(delta_t) + log(body_mass_g), V = vi, random = list(~1|study_ID, ~1|error), data = data_verts)
summary(model1)

# Fit model with phylogeny-ROM
model2_phylo <- rma.mv(yi ~ scale(delta_t) + log(body_mass_g), V = vi, random = list(~1|study_ID, ~1|species_rotl,  ~1|obs), R = list(species_rotl = PhyloA), data = data_verts_ROM)



# Fit model with phylogeny- CVR
model1_CVR <- rma.mv(yi ~ scale(delta_t) + log(body_mass_g), V = vi, random = list(~1|study_ID, ~1|species_rotl,  ~1|obs), R = list(species_rotl = PhyloA), data = data_verts_CVR)
























# Arm-based
arms_data <- read.csv("./data/lab_dive_durations_armbased.csv")

arms_data2 <- arms_data %>%
			filter(!study_ID == 13) %>%
			filter(!is.na(body_mass_g))

    arms_data2$mean_sv <- with(arms_data2, m_sv(mean, sd, n))
	  arms_data2$sd_sv <- with(arms_data2, sd_sv(n))

prior_slope <- list(R = list(V = 1, nu = 0.001),
                    G = list(G1 = list(V=1, nu = 0.02),
                             G2 = list(V = diag(2), nu = 2)))


A <- inverseA(tree, nodes = "TIPS", scale = TRUE)$Ainv
arms_data2$species_rotl <- paste0(arms_data2$genus, "_", arms_data2$species_new)

# Check what is different
setdiff(unique(arms_data2$species_rotl), sort(rownames(PhyloA)))

# Fix what is different in data
data_verts$species_rotl <- ifelse(data_verts$species_rotl == "Chrysemys_dorbignyi", "Trachemys_dorbigni", data_verts$species_rotl)
data_verts$species_rotl <- ifelse(data_verts$species_rotl == "Triturus_alpestris", "Ichthyosaura_alpestris", data_verts$species_rotl)





model6 <- MCMCglmm(log(mean) ~ log(body_mass_g) + T, mev = arms_data2$mean_sv, random = ~us(1):study_ID + us(1+T):species, data = arms_data2, prior = prior_slope, nitt = 50000, burnin = 10000, thin = 30)
summary(model6)


model7 <- MCMCglmm(log(sd) ~ log(mean) + log(body_mass_g) + T, mev = arms_data2$sd_sv, random = ~us(1):study_ID + us(1+T):species_rotl, ginverse = list(species_rotl = A), data = arms_data2, prior = prior_slope, nitt = 50000, burnin = 10000, thin = 30)
summary(model7)

ginverse = list(species_rotl = A)



# corrrelation = cov(X,Y) / (sd(X)*sd(Y))

cov_int_slope <- model6$VCV[,"T:(Intercept).species"]
       var_in <- model6$VCV[,"(Intercept):(Intercept).species"]
    var_slope <- model6$VCV[,"T:T.species"]

r <- cov_int_slope / (sqrt(var_in)*sqrt(var_slope))
mean(r)
HPDinterval(r)
hist(r)
quantile(r, probs =  c(0.25, 0.75))


