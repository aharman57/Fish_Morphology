library(tidyverse)
library(corrplot)
library(broom.mixed)
library(dotwhisker)
library(ggplot2); theme_set(theme_bw())
library(car)
library(geomorph)
library(effects)

Morph <- read_csv("Morph_Data_2016-2017.csv")
Morph_clean_body <- (Morph
                %>% dplyr::select(-c(1:2, 4, 9:12, 15:17, 20:25))
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
Morph_clean_body$Treatment <- as.factor(Morph_clean_body$Treatment)
Morph_clean_body$age <- as.factor(Morph_clean_body$age)

save(Morph_clean_body, file = "Morph_clean.R")

par(mfrow=c(1,1),mar=c(2,3,1.5,1),mgp=c(2,1,0))

#covariance matrix of response traits:
cov(Morph_clean_body[,1:7])
cormatrix <- cor(Morph_clean_body[,1:7])
corrplot.mixed(cormatrix, lower="ellipse", upper="number")

pairs(Morph_clean_body[, 1:7],
      pch = ".", gap = 0)

#logging variables:
Morph_log <- (Morph_clean_body
              %>% mutate(Length = log(Length),
                         Eye = log(Eye),
                         Fin_Anterior = log(Fin_Anterior),
                         Fin_Min = log(Fin_Min),
                         Fin_Posterior = log(Fin_Posterior),
                         Jaw = log(Jaw),
                         Body_Weight = log(Body_Weight)
              )
)

save(Morph_log, file = "Morph_log.R")

pairs(Morph_log[, 1:7],
      pch = ".", gap = 0)

## New Plots ##
scatterplotMatrix( ~ Length + Eye + Jaw + Fin_Min + Fin_Anterior + Fin_Posterior + Body_Weight,
                   ellipse = list(fill=TRUE, fill.alpha=0.6), data = Morph_log, gap = 0, regLine=FALSE, smooth=FALSE,
                   plot.points = F, pch = 20, cex  = 0.5, col=c("grey30", "grey0", "grey80"), groups=Morph_log$Treatment, by.groups=TRUE)

#eigenvalues after logging:
eig_vals_log <- svd(cov(Morph_log[, 1:7]))$d
prod(eig_vals_log)

#diagnostic plots for separate models:
model_length <- lm(Length ~ Treatment*age, data=Morph_log)
par(mfrow=c(2,2),mar=c(2,3,1.5,1),mgp=c(2,1,0))
plot(model_length, main = "Log Transformed Length")
model_eye <- lm(Eye ~ Treatment*age, data=Morph_log)
plot(model_eye, main = "Log Transformed Eye")
model_fin_anterior <- lm(Fin_Anterior ~ Treatment*age, data=Morph_log)
plot(model_fin_anterior, main = "Log Transformed Fin Anterior")
model_fin_min <- lm(Fin_Min ~ Treatment*age, data=Morph_log)
plot(model_fin_min, main = "Log Transformed Fin Min")
model_fin_posterior <- lm(Fin_Posterior ~ Treatment*age, data = Morph_log)
plot(model_fin_posterior, main = "Log Transformed Fin Posterior")
model_jaw <- lm(Jaw ~ Treatment*age, data=Morph_log)
plot(model_jaw, main = "Log Transformed Jaw")
model_weight <- lm(Body_Weight ~ Treatment*age, data = Morph_log)
plot(model_weight, main = "Log Transformed Body Weight")

library(MASS) #putting this here because if we put it at the top, it masks select function from tidyverse and messes up data cleaning

#trying box-cox transformation instead of logging (using model as object):
BoxLength <- boxcox(Length ~ Treatment*age, data=Morph_clean_body, lambda = seq(-4,4,0.1))
Cox_Length <- data.frame(BoxLength$x, BoxLength$y)
Cox_Length <- arrange(Cox_Length, desc(BoxLength.y))
lambda_length = Cox_Length[1, "BoxLength.x"]

BoxEye <- boxcox(Eye ~ Treatment*age, data=Morph_clean_body, lambda = seq(-4,4,0.1))
Cox_Eye = data.frame(BoxEye$x, BoxEye$y)
Cox_Eye <- arrange(Cox_Eye, desc(BoxEye.y))
lambda_eye = Cox_Eye[1, "BoxEye.x"]

BoxFinAnt <- boxcox(Fin_Anterior ~ Treatment*age, data=Morph_clean_body, lambda = seq(-4,4,0.1))
Cox_FinAnt = data.frame(BoxFinAnt$x, BoxFinAnt$y)
Cox_FinAnt <- arrange(Cox_FinAnt, desc(BoxFinAnt.y))
lambda_FinAnt = Cox_FinAnt[1, "BoxFinAnt.x"]

BoxFinMin <- boxcox(Fin_Min ~ Treatment*age, data=Morph_clean_body, lambda = seq(-4,4,0.1))
Cox_FinMin = data.frame(BoxFinMin$x, BoxFinMin$y)
Cox_FinMin <- arrange(Cox_FinMin, desc(BoxFinMin.y))
lambda_FinMin = Cox_FinMin[1, "BoxFinMin.x"]

BoxFinPost <- boxcox(Fin_Posterior ~ Treatment*age, data=Morph_clean_body, lambda = seq(-4,4,0.1))
Cox_FinPost = data.frame(BoxFinPost$x, BoxFinPost$y)
Cox_FinPost <- arrange(Cox_FinPost, desc(BoxFinPost.y))
lambda_FinPost = Cox_FinPost[1, "BoxFinPost.x"]

BoxJaw <- boxcox(Jaw ~ Treatment*age, data=Morph_clean_body, lambda = seq(-4, 4, 0.1))
Cox_Jaw = data.frame(BoxJaw$x, BoxJaw$y)
Cox_Jaw <- arrange(Cox_Jaw, desc(BoxJaw.y))
lambda_Jaw = Cox_Jaw[1, "BoxJaw.x"]

BoxWeight <- boxcox(Body_Weight ~ Treatment*age, data=Morph_clean_body, lambda = seq(-4, 4, 0.1))
Cox_Weight = data.frame(BoxWeight$x, BoxWeight$y)
Cox_Weight <- arrange(Cox_Weight, desc(BoxWeight.y))
lambda_Weight = Cox_Weight[1, "BoxWeight.x"]

Morph_Box <- (Morph_clean_body
              %>% mutate(Length = (Length ^ lambda_length - 1)/lambda_length)
              %>% mutate(Eye = (Eye ^ lambda_eye - 1)/lambda_eye)
              %>% mutate(Fin_Anterior = (Fin_Anterior ^ lambda_FinAnt - 1)/lambda_FinAnt)
              %>% mutate(Fin_Min = (Fin_Min ^ lambda_FinMin - 1)/lambda_FinMin)
              %>% mutate(Fin_Posterior = (Fin_Posterior ^ lambda_FinPost - 1)/lambda_FinPost)
              %>% mutate(Jaw = (Jaw ^ lambda_Jaw - 1)/lambda_Jaw)
              %>% mutate(Body_Weight = (Body_Weight ^ lambda_Weight - 1)/lambda_Weight)
)

#models and diagnostics with Box-Cox transformation:
par(mfrow=c(2,2),mar=c(2,3,1.5,1),mgp=c(2,1,0))
model_lengthbox <- lm(Length ~ Treatment*age, data=Morph_Box)
plot(model_lengthbox, main = "Box Transformed Length")
model_eyebox <- lm(Eye ~ Treatment*age, data=Morph_Box)
plot(model_eyebox, main = "Box Transformed Eye")
model_FinAntBox <- lm(Fin_Anterior ~ Treatment*age, data=Morph_Box)
plot(model_FinAntBox, main = "Box Transformed Fin Anterior")
model_FinMinBox <- lm(Fin_Min ~ Treatment*age, data=Morph_Box)
plot(model_FinMinBox, main = "Box Transformed Fin Min")
model_FinPostBox <- lm(Fin_Posterior ~ Treatment*age, data=Morph_Box)
plot(model_FinPostBox, main = "Box Transformed Fin Posterior")
model_jawbox <- lm(Jaw ~ Treatment*age, data=Morph_Box)
plot(model_jawbox, main = "Box Transformed Jaw")
model_weightbox <- lm(Body_Weight ~ Treatment*age, data=Morph_Box)
plot(model_weightbox, main = "Box Transformed Body Weight")

#multivariate model
mlm_fit1_log <- lm(as.matrix(Morph_log[,1:7]) ~ Treatment*age, data = Morph_log)
summary(manova(mlm_fit1_log), test = "Wilks")
coef(mlm_fit1_log)
exp(coef(mlm_fit1_log)) #back-transform to get biologically relevant effects

#magnitude of treatment and age constrast vectors - but what does this really mean?
sqrt(t(coef(mlm_fit1_log)[2,]) %*% coef(mlm_fit1_log)[2,])
sqrt(t(coef(mlm_fit1_log)[3,]) %*% coef(mlm_fit1_log)[3,])
sqrt(t(coef(mlm_fit1_log)[4,]) %*% coef(mlm_fit1_log)[4,])

#code for coefficient of determination:
sum(diag(cov(mlm_fit1_log$fitted)))/sum(diag(cov(Morph_log[,1:7])))
#model accounts for 66% of variance

#visualization:
plot(allEffects(mlm_fit1_log)) #this sort of works - maybe try to fix it up a bit - issue with ages
#is a ggplot possible?

#geomorph model:
mlm_fit2_log <- procD.lm(f1 = Morph_log[, 1:7] ~ Treatment*age, 
                     data = Morph_log, iter = 5000 )
summary(mlm_fit2_log)
coef(mlm_fit2_log)
#this basically gives same answer as first model

#trying a permutation test like in Ian's paper:

body_treatment_perm <- rep( NA, 1999 )
for(i in 1:1999){ 
  body_treatment_perm[i] <- summary( manova(lm( as.matrix( Morph_log[ sample(nrow(Morph_log), nrow(Morph_log), replace=F) ,1:7] ) ~ Morph_log$Treatment*Morph_log$age ) ))$stats[1,2]}
hist(body_treatment_perm, xlim=c(-1,1))
abline( v=summary( manova( mlm_fit1_log ))$stats[1,2], col="red")
#pseudo-p-val
mean(c(body_treatment_perm >= summary( manova( mlm_fit1_log ))$stats[1,2], 1))

#same, testing age:
body_age_perm <- rep( NA, 1999 )
for(i in 1:1999){ 
  body_age_perm[i] <- summary( manova(lm( as.matrix( Morph_log[ sample(nrow(Morph_log), nrow(Morph_log), replace=F) ,1:7] ) ~ Morph_log$Treatment*Morph_log$age ) ))$stats[2,2]}
hist(body_age_perm, xlim=c(-1.5,1.5))
abline( v=summary( manova( mlm_fit1_log ))$stats[2,2], col="red")
#pseudo-p-val
mean(c(body_age_perm >= summary( manova( mlm_fit1_log ))$stats[2,2], 1))

#same, testing interaction:
body_interact_perm <- rep( NA, 1999 )
for(i in 1:1999){ 
  body_interact_perm[i] <- summary( manova(lm( as.matrix( Morph_log[ sample(nrow(Morph_log), nrow(Morph_log), replace=F) ,1:7] ) ~ Morph_log$Treatment*Morph_log$age ) ))$stats[3,2]}
hist(body_interact_perm, xlim=c(-0.5,0.5))
abline( v=summary( manova( mlm_fit1_log ))$stats[3,2], col="red")
#pseudo-p-val
mean(c(body_interact_perm >= summary( manova( mlm_fit1_log ))$stats[3,2], 1))
