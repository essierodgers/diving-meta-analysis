#load libraries
install.packages("metafor")
install.packages("metagear", dependencies=TRUE)
library(metagear)
library(metafor)
library(robumeta)
library(dplyr)
library(MAd)
library(clubSandwich)
library(MuMIn)
install.packages("pacman") # install this on new computers/updated R
pacman::p_load(knitr, readxl, metafor, rotl, ape, Matrix, dplyr, purrr, stringr)
install.packages("BiocManager") 
BiocManager::install("EBImage")
library(EBImage)
library(knitr)
library(readxl)
library(rotl)
library(stringr)
install.packages("installr")
library(installr)
updateR()

#Activity rate data:
attach(activity_all2)
activity_all2 <- escalc(m1i = m2i, sd1i = sd2i, n1i= n2i, m2i = m1i, sd2i = sd1i, n2i=n1i, 
                      data = activity_all2, measure = "ROM", 
                      append = TRUE)
activity
options(max.print=1000000)
activity

str(activity_all2)
activity_all2$study_id = as.character(activity_all2$study_id)
activity_all2 <- as.data.frame(activity_all2)


activity_all2$RR <- log(activity_all2$m2i/activity_all2$m1i)
round(activity_all2$RR, 3)
V = covariance_commonControl(activity_all2, "study_id", "m2i", "sd2i", "n2i", "m1i", "sd1i", "n1i", metric = "RR")
round(V[[1]], 3)
activity_all2

V <- bldiag(lapply(split(activity[,c("vi")], activity$study_id), as.matrix))


suppressWarnings(suppressMessages(library(metafor))) # remove all messages when loading package 
theCovarianceMatrix <- V[[1]] 
theAlignedData <- V[[2]]
V <- bldiag(V)

activity1 <- rma.mv(RR, # a simple model that only pools the 3 effect sizes
       V=theCovarianceMatrix,  # inclusion of the sample VCV matrix
       random = list( ~1 |study_id,~1|species_name), 
       data = activity_all2, # the dataset with the effect sizes
       method = "REML",
       digits = 3)
print(activity1, digits = 3)

activity3 <- rma(RR, # a simple model that only pools the 3 effect sizes
                    vi,  # inclusion of the sample VCV matrix
                    data = activity_all2, # the dataset with the effect sizes
                    method = "REML",
                    digits = 3)
qqnorm(activity1, main = "Mixed-Effects Model")
qqnorm(activity3)

qqnorm(resid(activity1))
qqline(resid(activity1))



#Publication bias:
funnel(activity1, level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0) #Forest plot: 1st test of publication bias
fsn(yi, vi, data=activity, type="Rosenberg") #Fail safe number: 3rd test of publication bias


#random effects model
activitym1 <- rma.mv(yi,vi, random = ~1|study_id, struct="UN", data=activity)
print(activitym1, digits = 3)
profile(activity1)

#significant effect of nitrate on the activity of aquatic taxa (-0.814, CI: -0.972, -0.65)
#significant heterogeneity (Q 42 = 81.9, P <  0.001***) meaning that there is significant variation in the effect size

(exp(-0.29)-1)*100
#An effect size of 0.814 (+/- 0.081) indicates a 125.7% (+/- 8.4%, se) decrease 
#in activity


#Publication bias:
funnel(activitym1, level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0) #Forest plot: 1st test of publication bias
regtest(activitym1)
trimfill(activitym1) #Trimfill method: test and adjusts for publication bias
funnel(trimfill(activitym3), level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0)
fsn(yi, vi, data=activity, type="Rosenberg") #Fail safe number: 3rd test of publication bias

#Radial plot- simialar tool to visualise the effect size variablility
radial(activity1, main = "Random effects model")


#Multivariate analysis_Activity_all
activity_multi1 <-rma.mv(RR,  V=theCovarianceMatrix, 
                        mods= ~ logconcentration + logtime + factor(nitrate) + factor(Group) +  Temperature, 
                        random = list(~1|species_name, ~1|study_id), 
                        digits = 4, 
                        data=activity_all2, 
                        method="REML")
summary(activity_multi1)

install.packages("glmulti")
install.packages("rJava")
library(rJava)
library(glmulti)


rma.glmulti <- function(formula, data, ...)
  rma.mv(formula, vi, data=data, method="ML", ...)
active1 <- glmulti(yi ~ logconcentration + logtime + factor(nitrate) +Temperature + factor(taxa) + factor(life.history) + factor(source), data=activity, level=1, fitfunction=rma.glmulti, crit="aicc")
print(active1)
plot(active1)
tmp <- weightable(active1)
tmp <- tmp[tmp$aicc <= min(tmp$aicc),]
tmp

summary(active1@objects[[1]])

#Full model
activity_multi1 <-rma.mv(RR,  V=theCovarianceMatrix, 
                         mods= ~ logconcentration + logtime + factor(nitrate) + factor(Group) +  factor(stage) + factor(source) + Temperature, 
                         random = list(~1|species_name, ~1|study_id), 
                         R = list(species_name = phylo_cor),
                         digits = 4, 
                         data=activity_all2, 
                         method="FE")
summary(activity_multi1)


#Drop source and stage
activity_multi2 <-rma.mv(RR,  V=theCovarianceMatrix, 
                         mods= ~ logconcentration + logtime + factor(nitrate) + factor(Group) + Temperature, 
                         random = list(~1|species_name, ~1|study_id), 
                         R = list(species_name = phylo_cor),
                         digits = 4, 
                         data=activity_all2, 
                         method="FE")
summary(activity_multi2)

#Drop Group
activity_multi3 <-rma.mv(RR,  V=theCovarianceMatrix, 
                         mods= ~ logconcentration + logtime + factor(nitrate) + Temperature, 
                         random = list(~1|species_name, ~1|study_id), 
                         R = list(species_name = phylo_cor),
                         digits = 4, 
                         data=activity_all2, 
                         method="FE")
summary(activity_multi3)

#Drop temperature 
activity_multi4 <-rma.mv(RR,  V=theCovarianceMatrix, 
                         mods= ~ logconcentration + logtime + factor(nitrate), 
                         random = list(~1|species_name, ~1|study_id), 
                         R = list(species_name = phylo_cor),
                         digits = 4, 
                         data=activity_all2, 
                         method="FE")
summary(activity_multi4)

#Drop nitrate 
activity_multi5 <-rma.mv(RR,  V=theCovarianceMatrix, 
                         mods= ~ logconcentration + logtime + Temperature, 
                         random = list(~1|species_name, ~1|study_id), 
                         R = list(species_name = phylo_cor),
                         digits = 4, 
                         data=activity_all2, 
                         method="FE")
summary(activity_multi5)


anova(activity_multi1, activity_multi2, activity_multi3, activity_multi4, activity_multi5)

AIC(activity_multi1,activity_multi2, activity_multi5)


#Using a Robust variance estimation (RVE) method which  accounts for dependent effect sizes.  




#publication bias:
funnel(activity_multi, level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0) #Forest plot: 1st test of publication bias
regtest(activity_multi, predictor="ni")
fsn(yi, vi, data=activity, type="Rosenberg") #Fail safe number: 3rd test of publication bias

activity_taxa <-rma.mv(RR, V=theCovarianceMatrix,
                    mods = ~Group-1 ,
                    random = list(~1|species_name, ~1|study_id), 
                    R = list(species_name = phylo_cor),
                    digits = 4, 
                    data=activity_all2, 
                    method="FE")
activity_taxa


activity_nitrate <-rma.mv(RR, V=theCovarianceMatrix,
                          mods = ~nitrate-1 ,
                          random = list(~1|species_name, ~1|study_id), 
                          R = list(species_name = phylo_cor),
                          digits = 4, 
                          data=activity_all2, 
                          method="FE")
activity_nitrate

activity_lifehistory <-rma.mv(RR, V=theCovarianceMatrix,
                              mods = ~stage-1 ,
                              random = list(~1|species_name, ~1|study_id), 
                              R = list(species_name = phylo_cor),
                              digits = 4, 
                              data=activity_all2, 
                              method="FE")
activity_lifehistory

activity_source <-rma.mv(RR, V=theCovarianceMatrix,
                         mods = ~source-1 ,
                         random = list(~1|species_name, ~1|study_id), 
                         R = list(species_name = phylo_cor),
                         digits = 4, 
                         data=activity_all2, 
                         method="FE")
activity_source

activity_concen <-rma.mv(RR, V=theCovarianceMatrix,
                         mods = ~logconcentration ,
                         random = list(~1|species_name, ~1|study_id), 
                         R = list(species_name = phylo_cor),
                         digits = 4, 
                         data=activity_all2, 
                         method="FE")
activity_concen

activity_time <-rma(yi, vi, mods = ~logtime-1 , data=activity, method="ML" )
activity_time

activity_temperature <-rma(yi, vi, mods = ~Temperature-1 , data=activity, method="ML" )
activity_temperature


#Publication bias:
funnel(activity_multi, level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0) #First test of publication bias
regtest(activity_multi, model="lm")
fsn(yi, vi, data=activity, type="Rosenberg") #Fail safe number: 3rd test of publication bias


install.packages("ggplot2")
library(ggplot2)


ggplot(activity, aes(x= logconcentration, y=yi, size=n1i))+
  geom_point(alpha=0.25) +
  geom_smooth(method=lm , color="black")+ 
  scale_size_continuous(range = c(2, 12))+
  xlab("Log concentration")+
  ylab("Effect size (lnRR)")+
  labs( size = "Sample size (n)" ) +
  theme_classic()+
  theme(
    text = element_text(size=12),
    legend.justification = c("left", "bottom"),
    panel.border = element_blank()
    )+
  geom_hline(yintercept=0, linetype="dashed", color="black", size=0.5)


ggplot(activity, aes(x= logtime, y=yi, size=n1i))+
  geom_point(alpha=0.25) +
  geom_smooth(method=lm , color="black")+ 
  scale_size_continuous(range = c(2, 12))+
  xlab("Log duration")+
  ylab("Effect size (lnRR)")+
  labs( size = "Sample size (n)" ) +
  theme_classic()+
  theme(
    text = element_text(size=12),
    legend.justification = c("left", "bottom"),
    panel.border = element_blank()
    )+
  geom_hline(yintercept=0, linetype="dashed", color="black", size=0.5)


model1 = lm(yi~logtime, data=activity)
summary(model1)

model1 = lm(yi~logconcentration, data=activity)
summary(model1)
