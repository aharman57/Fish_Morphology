library(tidyverse)
library(corrplot)
library(broom.mixed)
library(ggplot2); theme_set(theme_bw())
library(car)
library(geomorph)
library(effects)

Morph <- read_csv("Morph_Data_2016-2017.csv")
Morph_clean_No_8C <- (Morph
                     %>% dplyr::select(-c(1:2, 4, 9:12, 15:17, 20:25)) #got rid of other variables we probably won't use
                     %>% rename(Length = "Length (mm)",
                                Eye = "Eye size (mm)",
                                Fin_Anterior = "Measure Dorsal Fin Anterior Maximum",
                                Fin_Min = "Measure Dorsal Fin Minimum",
                                Fin_Posterior = "Measure Dorsal Fin Posterior Maximum",
                                Jaw = "Jaw gape (um)",
                                Body_Weight = "body weight",
                                Treatment = "treatment group"
                     )
                     %>% filter(Treatment != 8)
                     %>% na.omit()
)

Morph_clean_No_8C$Fin_Anterior <- as.numeric(Morph_clean_No_8C$Fin_Anterior)
Morph_clean_No_8C$Fin_Min <- as.numeric(Morph_clean_No_8C$Fin_Min)
Morph_clean_No_8C$Fin_Posterior <- as.numeric(Morph_clean_No_8C$Fin_Posterior)
Morph_clean_No_8C$Treatment <- as.factor(Morph_clean_No_8C$Treatment)
Morph_clean_No_8C$age <- as.factor(Morph_clean_No_8C$age)

save(Morph_clean_No_8C, file = "Morph_clean_No_8C.R")

par(mfrow=c(1,1),mar=c(2,3,1.5,1),mgp=c(2,1,0))

#covariance matrix of response traits:
cov(Morph_clean_No_8C[,1:7])
cormatrix <- cor(Morph_clean_No_8C[,1:7])
corrplot.mixed(cormatrix, lower="ellipse", upper="number")

pairs(Morph_clean_No_8C[, 1:7],
      pch = ".", gap = 0)

#logging variables instead of scaling:
Morph_log_No_8C <- (Morph_clean_No_8C
              %>% mutate(Length = log(Length),
                         Eye = log(Eye),
                         Fin_Anterior = log(Fin_Anterior),
                         Fin_Min = log(Fin_Min),
                         Fin_Posterior = log(Fin_Posterior),
                         Jaw = log(Jaw),
                         Body_Weight = log(Body_Weight)
              )
)

save(Morph_log_No_8C, file = "Morph_log_No_8C.R")

pairs(Morph_log_No_8C[, 1:7],
      pch = ".", gap = 0)

## New Plots ##
scatterplotMatrix( ~ Length + Eye + Jaw + Fin_Min + Fin_Anterior + Fin_Posterior + Body_Weight,
                   ellipse = list(fill=TRUE, fill.alpha=0.6), data = Morph_log_No_8C, gap = 0, regLine=FALSE, smooth=FALSE,
                   plot.points = F, pch = 20, cex  = 0.5, col=c("grey30", "grey0", "grey80"), groups=Morph_log_No_8C$Treatment, by.groups=TRUE)


#diagnostic plots for separate models:
model_length <- lm(Length ~ Treatment*age, data=Morph_log_No_8C)
par(mfrow=c(2,2),mar=c(2,3,1.5,1),mgp=c(2,1,0))
plot(model_length, main = "Log Transformed Length")
model_eye <- lm(Eye ~ Treatment*age, data=Morph_log_No_8C)
plot(model_eye, main = "Log Transformed Eye")
model_fin_anterior <- lm(Fin_Anterior ~ Treatment*age, data=Morph_log_No_8C)
plot(model_fin_anterior, main = "Log Transformed Fin Anterior")
model_fin_min <- lm(Fin_Min ~ Treatment*age, data=Morph_log_No_8C)
plot(model_fin_min, main = "Log Transformed Fin Min")
model_fin_posterior <- lm(Fin_Posterior ~ Treatment*age, data = Morph_log_No_8C)
plot(model_fin_posterior, main = "Log Transformed Fin Posterior")
model_jaw <- lm(Jaw ~ Treatment*age, data=Morph_log_No_8C)
plot(model_jaw, main = "Log Transformed Jaw")
model_weight <- lm(Body_Weight ~ Treatment*age, data = Morph_log_No_8C)
plot(model_weight, main = "Log Transformed Body Weight")

## NEW/Better diagnostic plots
plot_model(model_length, type="diag", terms=c("age","Treatment"))
plot_model(model_eye, type="diag", terms=c("age","Treatment"))
plot_model(model_fin_anterior, type="diag", terms=c("age","Treatment"))
plot_model(model_fin_posterior, type="diag", terms=c("age","Treatment"))
plot_model(model_fin_min, type="diag", terms=c("age","Treatment"))
plot_model(model_jaw, type="diag", terms=c("age","Treatment"))
plot_model(model_weight, type="diag", terms=c("age","Treatment"))

#multivariate model
mlm_fit1_log_No_8C <- lm(as.matrix(Morph_log_No_8C[,1:7]) ~ Treatment*age, data = Morph_log_No_8C)
summary(manova(mlm_fit1_log_No_8C), test = "Wilks")
coef(mlm_fit1_log_No_8C)
exp(coef(mlm_fit1_log_No_8C)) #back-transform to get biologically relevant effects

#magnitude of treatment and age constrast vectors - but what does this really mean?
sqrt(t(coef(mlm_fit1_log_No_8C)[2,]) %*% coef(mlm_fit1_log_No_8C)[2,])
sqrt(t(coef(mlm_fit1_log_No_8C)[3,]) %*% coef(mlm_fit1_log_No_8C)[3,])
sqrt(t(coef(mlm_fit1_log_No_8C)[4,]) %*% coef(mlm_fit1_log_No_8C)[4,])

#code for coefficient of determination:
sum(diag(cov(mlm_fit1_log_No_8C$fitted)))/sum(diag(cov(Morph_log_No_8C[,1:7])))
#model accounts for 68% of variance? seems high

par(mfrow=c(1,1),mar=c(2,3,1.5,1),mgp=c(2,1,0))

#visualization:

plot(allEffects(mlm_fit1_log_No_8C)) 

#geomorph model (to validate coefficients):
LogCov <- cov(Morph_log_No_8C[,1:7])
mlm_fit2_log_No_8C <- procD.lm(f1 = Morph_log_No_8C[, 1:7] ~ Treatment*age, 
                         data = Morph_log_No_8C, iter = 5000 )
summary(mlm_fit2_log_No_8C)
coef(mlm_fit2_log_No_8C)

#trying a permutation test like in Ian's paper:

body_treatment_perm <- rep( NA, 1999 )
for(i in 1:1999){ 
  body_treatment_perm[i] <- summary( manova(lm( as.matrix( Morph_log_No_8C[ sample(nrow(Morph_log_No_8C), nrow(Morph_log_No_8C), replace=F) ,1:7] ) ~ Morph_log_No_8C$Treatment*Morph_log_No_8C$age ) ))$stats[1,2]}
hist(body_treatment_perm, xlim=c(-1,1))
abline( v=summary( manova( mlm_fit1_log_No_8C ))$stats[1,2], col="red")
#pseudo-p-val
mean(c(body_treatment_perm >= summary( manova( mlm_fit1_log_No_8C ))$stats[1,2], 1))

#same, testing age:
body_age_perm <- rep( NA, 1999 )
for(i in 1:1999){ 
  body_age_perm[i] <- summary( manova(lm( as.matrix( Morph_log_No_8C[ sample(nrow(Morph_log_No_8C), nrow(Morph_log_No_8C), replace=F) ,1:7] ) ~ Morph_log_No_8C$Treatment*Morph_log_No_8C$age ) ))$stats[2,2]}
hist(body_age_perm, xlim=c(-1.5,1.5))
abline( v=summary( manova( mlm_fit1_log_No_8C ))$stats[2,2], col="red")
#pseudo-p-val
mean(c(body_age_perm >= summary( manova( mlm_fit1_log_No_8C ))$stats[2,2], 1))

#same, testing interaction:
body_interact_perm <- rep( NA, 1999 )
for(i in 1:1999){ 
  body_interact_perm[i] <- summary( manova(lm( as.matrix( Morph_log_No_8C[ sample(nrow(Morph_log_No_8C), nrow(Morph_log_No_8C), replace=F) ,1:7] ) ~ Morph_log_No_8C$Treatment*Morph_log_No_8C$age ) ))$stats[3,2]}
hist(body_interact_perm, xlim=c(-0.5,0.5))
abline( v=summary( manova( mlm_fit1_log_No_8C ))$stats[3,2], col="red")
#pseudo-p-val
mean(c(body_interact_perm >= summary( manova( mlm_fit1_log_No_8C ))$stats[3,2], 1))
