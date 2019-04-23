#combining yolk and body size into one model
library(tidyverse)
library(corrplot)
library(broom.mixed)
library(car)
library(geomorph)
library(effects)
library(sjPlot)
library(snakecase)
library(MASS)

Morph <- read_csv("Morph_Data_2016-2017.csv")
Morph_clean_body_yolk <- (Morph
                     %>% dplyr::select(-c(1:13, 16:17, 20:25))
                     %>% rename(Body_Weight = "body weight",
                                Yolk_Weight = "yolk weight",
                                Treatment = "treatment group"
                     )
                     %>% na.omit()
                     %>% filter(age != 28, age != 14)
                     %>% mutate(Treatment = as.factor(Treatment),
                                age = as.factor(age))
)

save(Morph_clean_body_yolk, file="Morph_bodyyolk")

## scale instead of logging - 0 values
Morph_scale_body_yolk <- (Morph_clean_body_yolk
                        %>% mutate(Body_Weight = scale(Body_Weight),
                                       Yolk_Weight = scale(Yolk_Weight)
                                   )
)

## Pairs plot
pairs(Morph_scale_body_yolk[, 1:2],
      pch = ".", gap = 0)

## New Plots ##
scatterplotMatrix( ~ Body_Weight + Yolk_Weight,
                   ellipse = list(fill=TRUE, fill.alpha=0.6), data = Morph_scale_body_yolk, gap = 0, regLine=FALSE, smooth=FALSE,
                   plot.points = F, pch = 20, cex  = 0.5, col=c("grey30", "grey0", "grey80"), groups=Morph_scale_body_yolk$Treatment, by.groups=TRUE, xlim=c(-3,3), ylim=c(-3,3))

## MLM model fit
mlm_fit1_scale_yolkbody <- lm(as.matrix(Morph_scale_body_yolk[,1:2]) ~ Treatment*age, data = Morph_scale_body_yolk)
summary(manova(mlm_fit1_scale_yolkbody), test = "Wilks")
coef(mlm_fit1_scale_yolkbody)

#back-transforming coefficients:
(coef(mlm_fit1_scale_yolkbody)[,1])*sd(Morph_clean_body_yolk$Body_Weight)
(coef(mlm_fit1_scale_yolkbody)[,2])*sd(Morph_clean_body_yolk$Yolk_Weight)

#magnitude of treatment and age constrast vectors
sqrt(t(coef(mlm_fit1_scale_yolkbody)[2,]) %*% coef(mlm_fit1_scale_yolkbody)[2,])
sqrt(t(coef(mlm_fit1_scale_yolkbody)[3,]) %*% coef(mlm_fit1_scale_yolkbody)[3,])
sqrt(t(coef(mlm_fit1_scale_yolkbody)[4,]) %*% coef(mlm_fit1_scale_yolkbody)[4,])

#code for coefficient of determination:
sum(diag(cov(mlm_fit1_scale_yolkbody$fitted)))/sum(diag(cov(Morph_scale_body_yolk[,1:2])))

## Diagnostics
par(mfrow=c(2,2),mar=c(2,3,1.5,1),mgp=c(2,1,0))
lm_Body_Weight <- lm(Body_Weight ~ Treatment*age, data=Morph_scale_body_yolk)
plot(lm_Body_Weight, main = "Body_Weight")
lm_Yolk_Weight <- lm(Yolk_Weight ~ Treatment*age, data=Morph_scale_body_yolk)
plot(lm_Yolk_Weight, main = "Yolk_Weight")

## Different/Better Diagnostic Tests
plot_model(lm_Body_Weight, type="diag", terms=c("age","Treatment"))
plot_model(lm_Yolk_Weight, type="diag", terms=c("age","Treatment"))
### lots of heteroscadicity in the yolk weight - to be expected??? 

## All-Effects plot
plot(allEffects(mlm_fit1_scale_yolkbody))

# SJ plot - predictions
plot_model(lm_Body_Weight, type="pred", terms=c("age","Treatment"))
plot_model(lm_Yolk_Weight, type="pred", terms=c("age","Treatment"))

##Permutations using geomorph:
mlm_fit2_scale_yolkbody <- procD.lm(f1 = Morph_scale_body_yolk[, 1:2] ~ Treatment*age, 
                                data = Morph_scale_body_yolk, iter = 5000 )
summary(mlm_fit2_scale_yolkbody)
coef(mlm_fit2_scale_yolkbody)

##### Permutations #####

## Treatment
yolkbody_treatment_perm <- rep( NA, 1000 )
for(i in 1:1000){ 
  yolkbody_treatment_perm[i] <- summary( manova(lm( as.matrix( Morph_scale_body_yolk[ sample(nrow(Morph_scale_body_yolk), nrow(Morph_scale_body_yolk), replace=F) ,1:2] ) ~ Morph_scale_body_yolk$Treatment*Morph_scale_body_yolk$age ) ))$stats[1,2]}
hist(yolkbody_treatment_perm, xlim=c(-1,1))
abline( v=summary( manova( mlm_fit1_scale_yolkbody ))$stats[1,2], col="red")
#pseudo-p-val
mean(c(yolkbody_treatment_perm >= summary( manova( mlm_fit1_scale_yolkbody))$stats[1,2], 1))

## Age
yolkbody_age_perm <- rep( NA, 1000 )
for(i in 1:1000){ 
  yolkbody_age_perm[i] <- summary( manova(lm( as.matrix( Morph_scale_body_yolk[ sample(nrow(Morph_scale_body_yolk), nrow(Morph_scale_body_yolk), replace=F) ,1:2] ) ~ Morph_scale_body_yolk$Treatment*Morph_scale_body_yolk$age ) ))$stats[2,2]}
hist(yolkbody_age_perm, xlim=c(-1.5,1.5))
abline( v=summary( manova( mlm_fit1_scale_yolkbody ))$stats[2,2], col="red")
#pseudo-p-val
mean(c(yolkbody_age_perm >= summary( manova( mlm_fit1_scale_yolkbody ))$stats[2,2], 1)) ## same value as for treatment??

## Interaction
yolkbody_interact_perm <- rep( NA, 1000 )
for(i in 1:1000){ 
  yolkbody_interact_perm[i] <- summary( manova(lm( as.matrix( Morph_scale_body_yolk[ sample(nrow(Morph_scale_body_yolk), nrow(Morph_scale_body_yolk), replace=F) ,1:2] ) ~ Morph_scale_body_yolk$Treatment*Morph_scale_body_yolk$age ) ))$stats[3,2]}
hist(yolkbody_interact_perm, xlim=c(-0.5,0.5))
abline( v=summary( manova( mlm_fit1_scale_yolkbody ))$stats[3,2], col="red")
#pseudo-p-val
mean(c(yolkbody_interact_perm >= summary( manova( mlm_fit1_scale_yolkbody ))$stats[3,2], 1)) ## same value again..