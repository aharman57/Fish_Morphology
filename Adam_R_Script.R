library(tidyverse)
library(MCMCglmm)
library(lme4)
library(corrplot)
library(broom.mixed)
library(dotwhisker)
library(ggplot2); theme_set(theme_bw())
library(car)
library(geomorph)
library(emmeans)
library(dotwhisker)
library(effects)

Morph <- read_csv("Morph_Data_2016-2017.csv")
Morph_clean_body <- (Morph
                     %>% select(-c(1:2, 4, 9:12, 15:17, 20:25)) #got rid of other variables we probably won't use
                     %>% rename(Length = "Length (mm)",
                                Eye = "Eye size (mm)",
                                Fin_Anterior = "Measure Dorsal Fin Anterior Maximum",
                                Fin_Min = "Measure Dorsal Fin Minimum",
                                Fin_Posterior = "Measure Dorsal Fin Posterior Maximum",
                                Jaw = "Jaw gape (um)",
                                Body_Weight = "body weight",
                                Treatment = "treatment group"
                     )
                     %>% na.omit()
)

Morph_clean_body$Fin_Anterior <- as.numeric(Morph_clean_body$Fin_Anterior)
Morph_clean_body$Fin_Min <- as.numeric(Morph_clean_body$Fin_Min)
Morph_clean_body$Fin_Posterior <- as.numeric(Morph_clean_body$Fin_Posterior)
Morph_clean_body$Treatment <- as.factor(Morph_clean_body$Treatment) ## changed to factor
#### maybe make the age a factor as well??
#Morph_clean_body$age <- as.factor(Morph_clean_body$age)

#covariance matrix of response traits:
cov(Morph_clean_body[,1:7])
cormatrix <- cor(Morph_clean_body[,1:7])
corrplot(cormatrix, method = "circle")

pairs(Morph_clean_body[, 1:7],
      pch = ".", gap = 0)

#logging variables instead of scaling:
Morph_log <- (Morph_clean_body
              %>% mutate(Length = log(Length),
                         Eye = log(Eye),
                         Fin_Anterior = log(Fin_Anterior),
                         Fin_Min = log(Fin_Min),
                         Fin_Posterior = log(Fin_Posterior),
                         Jaw = log(Jaw),
                         Body_Weight = log(Body_Weight),
              )
)

pairs(Morph_log[, 1:7],
      pch = ".", gap = 0)

#eigenvalues after logging:
eig_vals_log <- svd(cov(Morph_log[, 1:7]))$d
prod(eig_vals_log)
sum(eig_vals_log)

mlm_fit1_log <- lm(as.matrix(Morph_log[,1:7]) ~ Treatment*age, data = Morph_log)
plot(mlm_fit1_log)
summary(manova(mlm_fit1_log), test = "Wilks")
coef(mlm_fit1_log)
#would need to back-transform effect sizes to get to biologically relevant scale - see class lecture slides


###############################
######## Diagnostics ##########
###############################

#Length
lm_fit_log <- lm(Length ~ Treatment*age, data = Morph_log)
par(mfrow=c(2,2),mar=c(2,3,1.5,1),mgp=c(2,1,0))
plot(lm_fit_log, which=1:4)

#Body Weight
lm_fit_log <- lm(Body_Weight ~ Treatment*age, data = Morph_log)
par(mfrow=c(2,2),mar=c(2,3,1.5,1),mgp=c(2,1,0))
plot(lm_fit_log, which=1:4)

#Eye
lm_fit_log <- lm(Eye ~ Treatment*age, data = Morph_log)
par(mfrow=c(2,2),mar=c(2,3,1.5,1),mgp=c(2,1,0))
plot(lm_fit_log, which=1:4)

#Jaw - resid vs fitted is non-linear
lm_fit_log <- lm(Jaw ~ Treatment*age, data = Morph_log)
par(mfrow=c(2,2),mar=c(2,3,1.5,1),mgp=c(2,1,0))
plot(lm_fit_log, which=1:4)

#Fin_Anterior
lm_fit_log <- lm(Fin_Anterior ~ Treatment*age, data = Morph_log)
par(mfrow=c(2,2),mar=c(2,3,1.5,1),mgp=c(2,1,0))
plot(lm_fit_log, which=1:4)

#Fin_Posterior
lm_fit_log <- lm(Fin_Posterior ~ Treatment*age, data = Morph_log)
par(mfrow=c(2,2),mar=c(2,3,1.5,1),mgp=c(2,1,0))
plot(lm_fit_log, which=1:4)

#Fin_Min
lm_fit_log <- lm(Fin_Min ~ Treatment*age, data = Morph_log)
par(mfrow=c(2,2),mar=c(2,3,1.5,1),mgp=c(2,1,0))
plot(lm_fit_log, which=1:4)

#####_______________________________________#####

#magnitude of treatment and age constrast vectors - but what does this really mean?
sqrt(t(coef(mlm_fit1_log)[2,]) %*% coef(mlm_fit1_log)[2,])
sqrt(t(coef(mlm_fit1_log)[3,]) %*% coef(mlm_fit1_log)[3,])
sqrt(t(coef(mlm_fit1_log)[4,]) %*% coef(mlm_fit1_log)[4,])

#code for coefficient of determination:
sum(diag(cov(Morph_log[,1:7])))
sum(diag(cov(mlm_fit1_log$fitted)))
sum(diag(cov(mlm_fit1_log$fitted)))/sum(diag(cov(Morph_log[,1:7])))
#model accounts for 66% of variance? seems high

dwplot(mlm_fit1_log) #this one doesn't work for some reason
plot(allEffects(mlm_fit1_log)) #AFFECTED BY MAKING TREATMENT A FACTOR - looks different now

#is a ggplot possible?

#geomorph model:
mlm_fit2_log <- procD.lm(f1 = Morph_log[, 1:7] ~ Treatment*age, 
                         data = Morph_log, iter = 2000 )
summary(mlm_fit2_log)
coef(mlm_fit2_log)
#this basically gives same answer as first model

#create coefficient plots?
