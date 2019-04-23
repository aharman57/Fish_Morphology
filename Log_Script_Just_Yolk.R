library(tidyverse)
library(corrplot)
library(broom.mixed)
library(ggplot2); theme_set(theme_bw())
library(car)
library(geomorph)
library(effects)
library(MASS)

Morph <- read_csv("Morph_Data_2016-2017.csv")
Morph_clean_yolk <- (Morph
                     %>% dplyr::select(-c(1:9, 12:14, 16:17, 20:25))
                     %>% rename(Yolk_Width ="Measure Yolk Width",
                                Yolk_Height = "Measure Yolk Height",
                                Yolk_Weight = "yolk weight",
                                Treatment = "treatment group"
                     )
                     %>% mutate(Yolk_Width = (Yolk_Width/1000),
                                Yolk_Height = (Yolk_Height/1000),
                                Treatment = as.factor(Treatment),
                                age = as.factor(age)
                     )
                     %>% na.omit()
                     %>% filter(age != 28, age != 14)
)



save(Morph_clean_yolk, file = "Morph_yolk.R")

#covariance matrix of response traits:
cov(Morph_clean_yolk[,1:3])
cormatrix_yolk <- cor(Morph_clean_yolk[,1:3])
corrplot(cormatrix_yolk, method = "circle")

pairs(Morph_clean_yolk[, 1:3],
      pch = ".", gap = 0)

#scaling variables instead of logging - can't log because of zero values:
Morph_scale_yolk <- (Morph_clean_yolk
              %>% mutate(Yolk_Width = scale(Yolk_Width),
                         Yolk_Height = scale(Yolk_Height),
                         Yolk_Weight = scale(Yolk_Weight)
              )
)

pairs(Morph_scale_yolk[, 1:3],
      pch = ".", gap = 0)

#eigenvalues after logging - need to deal with zeros in weight first:
eig_vals_scale_yolk <- svd(cov(Morph_scale_yolk[, 1:3]))$d
prod(eig_vals_scale_yolk)

#diagnostic plots for separate linear models:
model_yolkwidth <- lm(Yolk_Width ~ Treatment*age, data= Morph_scale_yolk)
par(mfrow=c(2,2),mar=c(2,3,1.5,1),mgp=c(2,1,0))
plot(model_yolkwidth)
model_yolkheight <- lm(Yolk_Height ~ Treatment*age, data=Morph_scale_yolk)
plot(model_yolkheight)
model_yolkweight <- lm(Yolk_Weight ~ Treatment*age, data = Morph_scale_yolk)
plot(model_yolkweight) #some heteroscedasticity here...will use permutations to calculate p values

#multivariate model:
mlm_fit1_scale_yolk <- lm(as.matrix(Morph_scale_yolk[,1:3]) ~ Treatment*age, data = Morph_scale_yolk)
summary(manova(mlm_fit1_scale_yolk), test = "Wilks")
coef(mlm_fit1_scale_yolk)

#back-transforming coefficients:
(coef(mlm_fit1_scale_yolk)[,1])*sd(Morph_clean_yolk$Yolk_Width)
(coef(mlm_fit1_scale_yolk)[,2])*sd(Morph_clean_yolk$Yolk_Height)
(coef(mlm_fit1_scale_yolk)[,3])*sd(Morph_clean_yolk$Yolk_Weight)

#magnitude of treatment and age constrast vectors
sqrt(t(coef(mlm_fit1_scale_yolk)[2,]) %*% coef(mlm_fit1_scale_yolk)[2,])
sqrt(t(coef(mlm_fit1_scale_yolk)[3,]) %*% coef(mlm_fit1_scale_yolk)[3,])
sqrt(t(coef(mlm_fit1_scale_yolk)[4,]) %*% coef(mlm_fit1_scale_yolk)[4,])

#code for coefficient of determination:
sum(diag(cov(mlm_fit1_scale_yolk$fitted)))/sum(diag(cov(Morph_scale_yolk[,1:3])))
#model accounts for 52% of variance

#visualization:
plot(allEffects(mlm_fit1_scale_yolk))
 
#geomorph model:
mlm_fit2_scale_yolk <- procD.lm(f1 = Morph_scale_yolk[, 1:3] ~ Treatment*age, 
                         data = Morph_scale_yolk, iter = 5000 )
summary(mlm_fit2_scale_yolk)
coef(mlm_fit2_scale_yolk)
#this basically gives same answer as first model

##### Permutations #####

## Treatment
yolk_treatment_perm <- rep( NA, 1000 )
for(i in 1:1000){ 
  yolk_treatment_perm[i] <- summary( manova(lm( as.matrix( Morph_scale_yolk[ sample(nrow(Morph_scale_yolk), nrow(Morph_scale_yolk), replace=F) ,1:2] ) ~ Morph_scale_yolk$Treatment*Morph_scale_yolk$age ) ))$stats[1,2]}
hist(yolk_treatment_perm, xlim=c(-1,1))
abline( v=summary( manova( mlm_fit1_scale_yolk ))$stats[1,2], col="red")
#pseudo-p-val
mean(c(yolk_treatment_perm >= summary( manova( mlm_fit1_scale_yolk))$stats[1,2], 1))

## Age
yolk_age_perm <- rep( NA, 1000 )
for(i in 1:1000){ 
  yolk_age_perm[i] <- summary( manova(lm( as.matrix( Morph_scale_yolk[ sample(nrow(Morph_scale_yolk), nrow(Morph_scale_yolk), replace=F) ,1:2] ) ~ Morph_scale_yolk$Treatment*Morph_scale_yolk$age ) ))$stats[2,2]}
hist(yolk_age_perm, xlim=c(-1.5,1.5))
abline( v=summary( manova( mlm_fit1_scale_yolk ))$stats[2,2], col="red")
#pseudo-p-val
mean(c(yolk_age_perm >= summary( manova( mlm_fit1_scale_yolk ))$stats[2,2], 1)) ## same value as for treatment??

## Interaction
yolk_interact_perm <- rep( NA, 1000 )
for(i in 1:1000){ 
  yolk_interact_perm[i] <- summary( manova(lm( as.matrix( Morph_scale_yolk[ sample(nrow(Morph_scale_yolk), nrow(Morph_scale_yolk), replace=F) ,1:2] ) ~ Morph_scale_yolk$Treatment*Morph_scale_yolk$age ) ))$stats[3,2]}
hist(yolk_interact_perm, xlim=c(-0.5,0.5))
abline( v=summary( manova( mlm_fit1_scale_yolk ))$stats[3,2], col="red")
#pseudo-p-val
mean(c(yolk_interact_perm >= summary( manova( mlm_fit1_scale_yolk ))$stats[3,2], 1)) ## same value again..
