##Log(Length as a predictor)

library(tidyverse)
library(corrplot)
library(broom.mixed)
library(dotwhisker)
library(ggplot2); theme_set(theme_bw())
library(car)
library(geomorph)
library(emmeans)
library(dotwhisker)
library(effects)

Morph <- read_csv("Morph_Data_2016-2017.csv") #could probably remove this too
load("Morph_clean.R")
load("Morph_log.R")

#diagnostic plots for separate models:
par(mfrow=c(2,2),mar=c(2,3,1.5,1),mgp=c(2,1,0))
allom_model_eye <- lm(Eye ~ Treatment*age*Length, data=Morph_log)
plot(allom_model_eye, main = "Log Transformed Eye")  ##need to fix titles
allom_model_fin_anterior <- lm(Fin_Anterior ~ Treatment*age*Length, data=Morph_log)
plot(allom_model_fin_anterior, main = "Log Transformed Fin Anterior")
allom_model_fin_min <- lm(Fin_Min ~ Treatment*age*Length, data=Morph_log)
plot(allom_model_fin_min, main = "Log Transformed Fin Min")
allom_model_fin_posterior <- lm(Fin_Posterior ~ Treatment*age*Length, data = Morph_log)
plot(allom_model_fin_posterior, main = "Log Transformed Fin Posterior")
allom_model_jaw <- lm(Jaw ~ Treatment*age*Length, data=Morph_log)
plot(allom_model_jaw, main = "Log Transformed Jaw")
allom_model_weight <- lm(Body_Weight ~ Treatment*age*Length, data = Morph_log)
plot(allom_model_weight, main = "Log Transformed Body Weight")

#maybe instead of also showing box cox, say that we will use permutations to calculate p value instead
##because if we use box cox then interpreting it might be hard?
##tried box cox on one of the worse ones (fin anterior) and it did not improve and we got negative lambda value

library(MASS)
par(mfrow=c(1,1),mar=c(2,3,1.5,1),mgp=c(2,1,0))
BoxFinAnt <- boxcox(Fin_Anterior ~ Treatment*age*Length, data=Morph_clean_body, lambda = seq(-4,4,0.1))
Cox_FinAnt = data.frame(BoxFinAnt$x, BoxFinAnt$y)
Cox_FinAnt <- arrange(Cox_FinAnt, desc(BoxFinAnt.y))
lambda_FinAnt = Cox_FinAnt[1, "BoxFinAnt.x"]

Morph_box_test <- (Morph_clean_body
                   %>% mutate(Fin_Anterior = (Fin_Anterior ^ lambda_FinAnt - 1)/lambda_FinAnt))

allom_box_test <- lm(Fin_Anterior ~ Treatment*age*Length, data=Morph_box_test)
par(mfrow=c(2,2),mar=c(2,3,1.5,1),mgp=c(2,1,0))
plot(allom_box_test)

#multivariate model:
mlm_allom_log <- lm(as.matrix(Morph_log[,2:7]) ~ Treatment*age*Length, data = Morph_log)
summary(manova(mlm_allom_log), test = "Wilks")
coef(mlm_allom_log)
exp(coef(mlm_allom_log)) ##do we need to back-transform since length is on the log scale too? how to interpret?

##coefficient of determination:
sum(diag(cov(mlm_allom_log$fitted)))/sum(diag(cov(Morph_log[,2:7])))
#model accounts for ~ 77% of variance

par(mfrow=c(1,1),mar=c(2,3,1.5,1),mgp=c(2,1,0))

plot(allEffects(mlm_allom_log)) ##this doesnt work too well because too squished
#is there a way to do one or two at a time? or just do ggplot

#permutation test using geomorph
LogCov <- cov(Morph_log[,2:7])
mlm_allom_geo <- procD.lm(f1 = Morph_log[, 2:7] ~ Treatment*age*Length, 
                         data = Morph_log, iter = 5000 ) ##something not working here - fixed by removing specification of covariance matrix - assuming default is accurate??
summary(mlm_allom_geo)
coef(mlm_allom_geo)

#permutation test like in Ian's paper:

allom_treatment_perm <- rep( NA, 1999 )
allom_age_perm <- rep( NA, 1999)
allom_length_perm <- rep(NA, 1999)
allom_treatage_perm <- rep( NA, 1999)
allom_treatlen_perm <- rep( NA, 1999)
allom_agelen_perm <- rep( NA, 1999)
allom_interact_perm <- rep( NA, 1999)

for(i in 1:1999){ 
  allom_treatment_perm[i] <- summary( manova(lm( as.matrix( Morph_log[ sample(nrow(Morph_log), nrow(Morph_log), replace=F) ,2:7] ) ~ Morph_log$Treatment*Morph_log$age*Morph_log$Length ) ))$stats[1,2]
  allom_age_perm[i] <- summary( manova(lm( as.matrix( Morph_log[ sample(nrow(Morph_log), nrow(Morph_log), replace=F) ,2:7] ) ~ Morph_log$Treatment*Morph_log$age*Morph_log$Length ) ))$stats[2,2]
  allom_length_perm[i] <- summary( manova(lm( as.matrix( Morph_log[ sample(nrow(Morph_log), nrow(Morph_log), replace=F) ,2:7] ) ~ Morph_log$Treatment*Morph_log$age*Morph_log$Length ) ))$stats[3,2]
  allom_treatage_perm[i] <- summary( manova(lm( as.matrix( Morph_log[ sample(nrow(Morph_log), nrow(Morph_log), replace=F) ,2:7] ) ~ Morph_log$Treatment*Morph_log$age*Morph_log$Length ) ))$stats[4,2]
  allom_treatlen_perm[i] <- summary( manova(lm( as.matrix( Morph_log[ sample(nrow(Morph_log), nrow(Morph_log), replace=F) ,2:7] ) ~ Morph_log$Treatment*Morph_log$age*Morph_log$Length ) ))$stats[5,2]
  allom_agelen_perm[i] <- summary( manova(lm( as.matrix( Morph_log[ sample(nrow(Morph_log), nrow(Morph_log), replace=F) ,2:7] ) ~ Morph_log$Treatment*Morph_log$age*Morph_log$Length ) ))$stats[6,2]
  allom_interact_perm[i] <- summary( manova(lm( as.matrix( Morph_log[ sample(nrow(Morph_log), nrow(Morph_log), replace=F) ,2:7] ) ~ Morph_log$Treatment*Morph_log$age*Morph_log$Length ) ))$stats[7,2]}
hist(allom_treatment_perm, xlim=c(-1,1))
abline( v=summary( manova( mlm_allom_log ))$stats[1,2], col="red")
#pseudo-p-val
mean(c(body_treatment_perm >= summary( manova( mlm_allom_log ))$stats[1,2], 1))

