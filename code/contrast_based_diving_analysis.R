# Load libraries
pacman::p_load(metafor, MCMCglmm, tidyverse, rotl, phytools, corrplot)
library(devtools)
install_github("daniel1noble/metaAidR"); library(metaAidR)

# Read in the data
data <- read.csv("~/diving-meta-analysis/data/lab_dive_durations_contrast_highlow.csv")
str(data)
data$body_mass_g <- as.numeric(data$body_mass_g)

#separate verts from inverts
data_verts <- data %>% filter(study_ID != 13)


#how many species per study
dat <- data %>%
  group_by(study_ID) %>%
  summarise(n = length(unique(species))) 

#importing phylogeny from timetree
tree <- read.tree("~/diving-meta-analysis/data/order_data/vert_phylogeny.NWK")
plot(tree)
PhyloA <- vcv(tree, corr = TRUE)



# Separating genus and species
genus <- as.data.frame(do.call("rbind", str_split(str_trim(data_verts$species, side = "both"), " ")))
names(genus) <- c("genus", "species_new")
data_verts <- cbind(data_verts, genus)


data_verts$species_rotl <- paste0(data_verts$genus, "_", data_verts$species_new)




## Plots
plot_func <- function(data, x, y){
  plot(log(data[,y]) ~log(data[,x]), pch = 16, ylab = "log(y)", xlab = "log(x)") 
}

plot_func(data, "sd_control", "mean_control")


#calculating effect sizes, measure="ROM" means lnRR es, append=true means will add es to data set. yi= es estimate, vi sampling variance.
data_verts <- escalc(m1i = mean_control, m2i = mean_treatment, n1i = n_control, n2i = n_treatment, sd1i = sd_control, sd2i = sd_treatment, append = TRUE, measure ="ROM", data = data_verts)

# Check directionality – control on numerator
(log(data$mean_control) - log(data$mean_treatment) ) == (log(data$mean_control/data$mean_treatment))

#error calc-residual variance at the observation level
data_verts$error <- 1:dim(data_verts)[1]

# Meta-analytic multivariate model - wihtout intercept
model1 <- rma.mv(yi ~ -1 + scale(delta_t) + log(body_mass_g), V = vi, random = list(~1|study_ID, ~1|error), data = data)
summary(model1)

# Fit model with phylogeny
model2_phylo <- rma.mv(yi ~ scale(delta_t) + log(body_mass_g), V = vi, random = list(~1|study_ID, ~1|species_rotl,  ~1|error), R = list(species_rotl = PhyloA), data = data_verts)


# Check same number of levels
length(rownames(PhyloA))
length(unique(data_verts$species_rotl))

# Check what is different
setdiff(unique(data_verts$species_rotl), sort(rownames(PhyloA)))

# Fix what is different in data
data_verts$species_rotl <- ifelse(data_verts$species_rotl == "Chrysemys_dorbignyi", "Trachemys_dorbigni", data_verts$species_rotl)
data_verts$species_rotl <- ifelse(data_verts$species_rotl == "Triturus_alpestris", "Ichthyosaura_alpestris", data_verts$species_rotl)


