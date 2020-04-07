# Load libraries
pacman::p_load(metafor, MCMCglmm, tidyverse, rotl, phytools, corrplot)
library(devtools)
install_github("daniel1noble/metaAidR"); library(metaAidR)

#data processing
rerun = FALSE
if(rerun = TRUE){
  data <- read.csv("./data/lab_dive_durations_contrast.csv")
  data$body_mass_g <- as.numeric(data$body_mass_g)
  data$t_magnitude <- as.factor(data$t_magnitude)
  data$study_ID <- as.factor(data$study_ID)
  #separate verts from inverts
  data_verts <- data %>% filter(study_ID != 13)
  # Separating genus and species
  genus <- as.data.frame(do.call("rbind", str_split(str_trim(data_verts$species, side = "both"), " ")))
  names(genus) <- c("genus", "species_new")
  data_verts <- cbind(data_verts, genus)
  data_verts$species_rotl <- paste0(data_verts$genus, "_", data_verts$species_new)
  data_verts <- data_verts[,-31]
  #importing phylogeny from timetree
  tree <- read.tree("./data/order_data/vert_phylogeny.NWK")
  plot(tree)
  PhyloA <- vcv(tree, corr = TRUE)
  # Fix what is different in data
  data_verts$species_rotl <- ifelse(data_verts$species_rotl == "Chrysemys_dorbignyi", "Trachemys_dorbigni", data_verts$species_rotl)
  data_verts$species_rotl <- ifelse(data_verts$species_rotl == "Triturus_alpestris", "Ichthyosaura_alpestris", data_verts$species_rotl)
  write.csv(data_verts, "processed_data_file.csv")
}else{
  datafile <- read.csv("./processed_data_file.csv")
}


## Plots
plot_func <- function(data, x, y){
  plot(log(data[,y]) ~log(data[,x]), pch = 16, ylab = "log(y)", xlab = "log(x)") 
}

plot_func(data, "sd_control", "mean_control")


#calculating effect sizes, measure="ROM" means lnRR es, append=true means will add es to data set. yi= es estimate, vi sampling variance.
data_verts_ROM <- escalc(m1i = mean_treatment, m2i = mean_control, n1i = n_treatment, n2i = n_control, sd1i = sd_treatment, sd2i = sd_treatment, append = TRUE, measure ="ROM", data = data_verts)
data_verts_CVR <- escalc(m1i = mean_treatment, m2i = mean_control, n1i = n_treatment, n2i = n_control, sd1i = sd_treatment, sd2i = sd_treatment, append = TRUE, measure ="CVR", data = data_verts)


# Check directionality – control on numerator
(log(data$mean_control) - log(data$mean_treatment) ) == (log(data$mean_control/data$mean_treatment))

#error calc-residual variance at the observation level
data_verts_ROM$obs <- 1:dim(data_verts_ROM)[1]
data_verts_CVR$obs <- 1:dim(data_verts_CVR)[1]


# Make shared control matrix

make_VCV_matrix <- function(data, V, cluster, obs, type=c("vcv", "cor"), rho=0.5){
  
  if (missing(data)) 
    stop("Must specify dataframe via 'data' argument.")
  if (missing(V)) 
    stop("Must specify name of the variance variable via 'V' argument.")
  if (missing(cluster)) 
    stop("Must specify name of the clustering variable via 'cluster' argument.")
  if (missing(obs)) 
    obs <- 1:length(V)   
  if (missing(type)) 
    type <- "vcv" 
  
  new_matrix <- matrix(0,nrow = dim(data)[1],ncol = dim(data)[1]) #make empty matrix of the same size as data length
  rownames(new_matrix) <- data[ ,obs]
  colnames(new_matrix) <- data[ ,obs]
  # find start and end coordinates for the subsets
  shared_coord <- which(data[ ,cluster] %in% data[duplicated(data[ ,cluster]), cluster]==TRUE)
  # matrix of combinations of coordinates for each experiment with shared control
  combinations <- do.call("rbind", tapply(shared_coord, data[shared_coord,cluster], function(x) t(utils::combn(x,2))))
  
  if(type == "vcv"){
    # calculate covariance values between  values at the positions in shared_list and place them on the matrix
    for (i in 1:dim(combinations)[1]){
      p1 <- combinations[i,1]
      p2 <- combinations[i,2]
      p1_p2_cov <- rho * sqrt(data[p1,V]) * sqrt(data[p2,V])
      new_matrix[p1,p2] <- p1_p2_cov
      new_matrix[p2,p1] <- p1_p2_cov
    }
    diag(new_matrix) <- data[ ,V]   #add the diagonal
  }
  
  if(type == "cor"){
    # calculate covariance values between  values at the positions in shared_list and place them on the matrix
    for (i in 1:dim(combinations)[1]){
      p1 <- combinations[i,1]
      p2 <- combinations[i,2]
      p1_p2_cov <- rho
      new_matrix[p1,p2] <- p1_p2_cov
      new_matrix[p2,p1] <- p1_p2_cov
    }
    diag(new_matrix) <- 1   #add the diagonal of 1
  }
  
  return(new_matrix)
}

#Shared control for ROM dataset
data_verts_ROM$sc_cluster <- interaction(data_verts_ROM$study_ID, data_verts_ROM$shared_control)
V <- make_VCV_matrix(data_verts_ROM, "vi", "sc_cluster", type = "vcv", rho = 0.5)
corrplot(as.matrix(V))
write.csv(V, file = "sc_matrix_rom.csv")


#shared control for CVR dataset
data_verts_CVR$sc_cluster <- interaction(data_verts_CVR$study_ID, data_verts_CVR$shared_control)
V2 <- make_VCV_matrix(data_verts_CVR, "vi", "sc_cluster", type = "vcv", rho = 0.5)
corrplot(as.matrix(V))
write.csv(V, file = "sc_matrix_cvr.csv")


#FINAL MODELS- with shared control accounted for

#Mean models
#Overall effect of temperature increase on dive duration
model1 <- rma.mv(yi = yi, V = V, random = list(~1|study_ID, ~1|species_rotl,  ~1|obs), R = list(species_rotl = PhyloA), data = data_verts_ROM)
summary(model1)


#Model with moderators
model2 <- rma.mv(yi = yi, V = V, mods = ~ mean_t + delta_t + log(body_mass_g) + respiration_mode, random = list(~1|study_ID, ~1|species_rotl,  ~1|obs), R = list(species_rotl = PhyloA), data = data_verts_ROM)
summary(model2)

#Effect sizes for bimodal veresus aerial breathers
model3 <- rma.mv(yi = yi, V = V, mods = ~respiration_mode-1, random = list(~1|study_ID, ~1|species_rotl,  ~1|obs), R = list(species_rotl = PhyloA), data = data_verts_ROM)
summary(model3)

#Effect sizes for different magnitudes of temperature increase (+3, +5-7, +8-9, +>10C)
model4 <- rma.mv(yi = yi, V = V, mods = ~t_magnitude-1, random = list(~1|study_ID, ~1|species_rotl,  ~1|obs), R = list(species_rotl = PhyloA), data = data_verts_ROM)
summary(model4)


#Order effect sizes
model5 <- rma.mv(yi = yi, V = V, mods = ~order-1, random = list(~1|study_ID, ~1|species_rotl,  ~1|obs), R = list(species_rotl = PhyloA), data = data_verts_ROM)
summary(model5)


#CVR models
#Overall effect of temperature increase on dive duration variability
model1 <- rma.mv(yi = yi, V = V2, random = list(~1|study_ID, ~1|species_rotl,  ~1|obs), R = list(species_rotl = PhyloA), data = data_verts_CVR)
summary(model1)

#Model with moderators
model2 <- rma.mv(yi = yi, V = V2, mods = ~ mean_t + delta_t + log(body_mass_g) + respiration_mode, random = list(~1|study_ID, ~1|species_rotl,  ~1|obs), R = list(species_rotl = PhyloA), data = data_verts_CVR)
summary(model2)

#Effect sizes for bimodal veresus aerial breathers
model3 <- rma.mv(yi = yi, V = V2, mods = ~respiration_mode-1, random = list(~1|study_ID, ~1|species_rotl,  ~1|obs), R = list(species_rotl = PhyloA), data = data_verts_CVR)
summary(model3)

#Effect sizes for different magnitudes of temperature increase (+3, +5-7, +8-9, +>10C)
model4 <- rma.mv(yi = yi, V = V2, mods = ~t_magnitude-1, random = list(~1|study_ID, ~1|species_rotl,  ~1|obs), R = list(species_rotl = PhyloA), data = data_verts_CVR)
summary(model4)


#Order effect sizes
model5 <- rma.mv(yi = yi, V = V2, mods = ~order-1, random = list(~1|study_ID, ~1|species_rotl,  ~1|obs), R = list(species_rotl = PhyloA), data = data_verts_CVR)
summary(model5)





#OLD CODE-IGNORE


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
data_verts <- data_verts[,-31]


# Check same number of levels
length(rownames(PhyloA))
length(unique(data_verts$species_rotl))

# Check what is different
setdiff(unique(data_verts$species_rotl), sort(rownames(PhyloA)))

# Fix what is different in data
data_verts$species_rotl <- ifelse(data_verts$species_rotl == "Chrysemys_dorbignyi", "Trachemys_dorbigni", data_verts$species_rotl)
data_verts$species_rotl <- ifelse(data_verts$species_rotl == "Triturus_alpestris", "Ichthyosaura_alpestris", data_verts$species_rotl)


# Fit model with phylogeny-ROM
model1_phylo <- rma.mv(yi ~ mean_t + delta_t + log(body_mass_g) + respiration_mode, V = vi, random = list(~1|study_ID, ~1|species_rotl,  ~1|obs), R = list(species_rotl = PhyloA), data = data_verts_ROM)




model1_test <- rma.mv(yi = yi, V = vi, mods = ~ mean_t + delta_t + log(body_mass_g) + respiration_mode, random = list(~1|study_ID, ~1|species_rotl,  ~1|obs), R = list(species_rotl = PhyloA), data = data_verts_ROM)
model1_test <- rma.mv(yi = yi, V = vi, mods = ~respiration_mode-1, random = list(~1|study_ID, ~1|species_rotl,  ~1|obs), R = list(species_rotl = PhyloA), data = data_verts_ROM)




# Fit model with phylogeny- CVR
model1_CVR <- rma.mv(yi ~   mean_t + delta_t + log(body_mass_g), V = vi, random = list(~1|study_ID, ~1|species_rotl,  ~1|obs), R = list(species_rotl = PhyloA), data = data_verts_CVR)
model1_CVR <- rma.mv(yi = yi, V = vi,  random = list(~1|study_ID, ~1|species_rotl,  ~1|obs), R = list(species_rotl = PhyloA), data = data_verts_CVR)









# Meta-analytic multivariate model - wihtout intercept
model1 <- rma.mv(yi ~ -1 + scale(delta_t) + log(body_mass_g), V = vi, random = list(~1|study_ID, ~1|error), data = data_verts)
summary(model1)








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


