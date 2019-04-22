#combining yolk and body size into one model
#one of the packages is interfering with tidyverse... cant select data?? --> works now with these packages
library(tidyverse)
library(corrplot)
library(car)
library(geomorph)
library(effects)

Morph <- read_csv("Morph_Data_2016-2017.csv")
Morph_clean_body_yolk <- (Morph
                     %>% select(-c(1:13, 16:17, 20:25))
                     %>% rename(Body_Weight = "body weight",
                                Yolk_Weight = "yolk weight",
                                Treatment = "treatment group"
                     )
                     %>% na.omit()
                     %>% filter(age != 28, age != 14)
                     %>% mutate(Treatment = as.factor(Treatment),
                                age = as.factor(age))
)

## scale instead of logging - 0 values
Morph_scale_body_yolk <- (Morph_clean_body_yolk
                        %>% mutate(Body_Weight = scale(Body_Weight),
                                       Yolk_Weight = scale(Yolk_Weight)
                                   )
)

#################


mlm_fit1_scale_yolkbody <- lm(as.matrix(Morph_scale_body_yolk[,1:2]) ~ Treatment*age, data = Morph_scale_body_yolk)
summary(manova(mlm_fit1_scale_yolkbody), test = "Wilks")
coef(mlm_fit1_scale_yolkbody)

## Diagnostics
par(mfrow=c(2,2),mar=c(2,3,1.5,1),mgp=c(2,1,0))
lm_Body_Weight <- lm(Body_Weight ~ Treatment*age, data=Morph_scale_body_yolk)
plot(lm_Body_Weight, main = "Body_Weight")
lm_Yolk_Weight <- lm(Yolk_Weight ~ Treatment*age, data=Morph_scale_body_yolk)
plot(lm_Yolk_Weight, main = "Yolk_Weight")
### lots of heteroscadicity in the yolk weight - to be expected??? 

##Permutations using geomorph:
mlm_fit2_scale_yolkbody <- procD.lm(f1 = Morph_scale_body_yolk[, 1:2] ~ Treatment*age, 
                                data = Morph_scale_body_yolk, iter = 5000 )
summary(mlm_fit2_scale_yolkbody)
coef(mlm_fit2_scale_yolkbody)

##### Permutations #####

library(MASS)

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

## All-Effects plot
Yolk_Effects <- allEffects(mlm_fit1_scale_yolkbody)
plot(allEffects(mlm_fit1_scale_yolkbody))
# Separates the predictors, can't seem to separate the response... always faceted??
plot(effect(mod=mlm_fit1_scale_yolkbody, term = "age", residuals=TRUE))
plot(effect(mod=mlm_fit1_scale_yolkbody, term = "Treatment"))

library(sjPlot)
library(snakecase)
sjp.int(mlm_fit1_scale_yolkbody, swap.pred = T)
plot_model(mlm_fit1_scale_yolkbody, type="pred", terms=c("age","Treatment"))

plot_model(lm_Body_Weight, type="pred", terms=c("age","Treatment"))
