library(tidyverse)
library(MCMCglmm)
library(lme4)
library(corrplot)
library(broom.mixed)
library(dotwhisker)
library(ggplot2); theme_set(theme_bw())
library(car)
library(geomorph)

Morph <- read_csv("Morph_Data_2016-2017.csv")
Morph_clean_yolk <- (Morph
                     %>% select(-c(1:9, 12:14, 16:17, 20:25)) #got rid of other variables we probably won't use
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

#covariance matrix of response traits:
cov(Morph_clean_yolk[,1:3])
cormatrix_yolk <- cor(Morph_clean_yolk[,1:3])
corrplot(cormatrix_yolk, method = "circle") #fix something here

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
sum(eig_vals_scale_yolk)

#diagnostic plots for separate linear models:
model_yolkwidth <- lm(Yolk_Width ~ Treatment*age, data= Morph_scale_yolk)
par(mfrow=c(2,2),mar=c(2,3,1.5,1),mgp=c(2,1,0))
plot(model_yolkwidth)
model_yolkheight <- lm(Yolk_Height ~ Treatment*age, data=Morph_scale_yolk)
plot(model_yolkheight) #issue here with invalid value - missing value?
model_yolkweight <- lm(Yolk_Weight ~ Treatment*age, data = Morph_scale_yolk)
plot(model_yolkweight) #some heteroskedasticity here...will use permutations to calculate p values

#multivariate model:
mlm_fit1_scale_yolk <- lm(as.matrix(Morph_scale_yolk[,1:3]) ~ Treatment*age, data = Morph_scale_yolk)
summary(manova(mlm_fit1_scale_yolk), test = "Wilks")
coef(mlm_fit1_scale_yolk)
#how to back-transform effect sizes when they are scaled?
#now that treatment is a factor, shouldn't it show the different levels (contrasts from baseline?)

#magnitude of treatment and age constrast vectors - but what does this really mean? DELETE?
sqrt(t(coef(mlm_fit1_scale_yolk)[2,]) %*% coef(mlm_fit1_scale_yolk)[2,])
sqrt(t(coef(mlm_fit1_scale_yolk)[3,]) %*% coef(mlm_fit1_scale_yolk)[3,])
sqrt(t(coef(mlm_fit1_scale_yolk)[4,]) %*% coef(mlm_fit1_scale_yolk)[4,])

#code for coefficient of determination:
sum(diag(cov(Morph_scale_yolk[,1:3])))
sum(diag(cov(mlm_fit1_scale_yolk$fitted)))
sum(diag(cov(mlm_fit1_scale_yolk$fitted)))/sum(diag(cov(Morph_scale_yolk[,1:3])))
#model accounts for 52% of variance

#figure out if we need to do permutation test stuff to assess whether data conform to assumptions
#visualization:
dwplot(mlm_fit1_scale_yolk) #this one doesn't work for some reason
plot(allEffects(mlm_fit1_scale_yolk)) #this sort of works - maybe try to fix it up a bit
plot(emmeans(mlm_fit1_scale_yolk, ~Treatment)) #this is useless
#is a ggplot possible? crate object with allEffects and then use ggplot?
 
#geomorph model:
mlm_fit2_scale_yolk <- procD.lm(f1 = Morph_scale_yolk[, 1:3] ~ Treatment*age, 
                         data = Morph_scale_yolk, iter = 5000 )
summary(mlm_fit2_scale_yolk)
coef(mlm_fit2_scale_yolk)
#this basically gives same answer as first model

#create coefficient plots?

##### Permutations ##### 

library(MASS) #load here because it masks tidyverse 'select' function

## Treatment
yolk_treatment_perm <- rep( NA, 1000 )
for(i in 1:1000){ 
  yolk_treatment_perm[i] <- summary( manova(lm( as.matrix( Morph_scale_yolk[ sample(nrow(Morph_scale_yolk), nrow(Morph_scale_yolk), replace=F) ,1:2] ) ~ Morph_scale_yolk$Treatment*Morph_scale_yolk$age ) ))$stats[1,2]}
hist(yolk_treatment_perm, xlim=c(-1,1))
abline( v=summary( manova( mlm_fit1_scale_yolk ))$stats[1,2], col="red")
#pseudo-p-val
mean(c(yolk_treatment_perm >= summary( manova( mlm_fit1_scale_yolk))$stats[1,2], 1))

## Age ##should change the name of the dataset potentially for each of these and for each script in case we want to run them all at once
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
