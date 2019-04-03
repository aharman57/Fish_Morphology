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
                                Yolk_Height = (Yolk_Height/1000)
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

#logging variables instead of scaling:
Morph_log_yolk <- (Morph_clean_yolk
              %>% mutate(Yolk_Width = log(Yolk_Width+1),
                         Yolk_Height = log(Yolk_Height+1),
                         Yolk_Weight = log(Yolk_Weight+1)
              )
)

pairs(Morph_log_yolk[, 1:3],
      pch = ".", gap = 0)

#eigenvalues after logging - need to deal with zeros in weight first:
eig_vals_log_yolk <- svd(cov(Morph_log_yolk[, 1:3]))$d
prod(eig_vals_log_yolk)
sum(eig_vals_log_yolk)

mlm_fit1_log_yolk <- lm(as.matrix(Morph_log_yolk[,1:3]) ~ Treatment*age, data = Morph_log_yolk)
summary(manova(mlm_fit1_log_yolk), test = "Wilks")
coef(mlm_fit1_log_yolk)
#would need to back-transform effect sizes to get to biologically relevant scale - see class lecture slides

#magnitude of treatment and age constrast vectors - but what does this really mean?
sqrt(t(coef(mlm_fit1_log_yolk)[2,]) %*% coef(mlm_fit1_log_yolk)[2,])
sqrt(t(coef(mlm_fit1_log_yolk)[3,]) %*% coef(mlm_fit1_log_yolk)[3,])
sqrt(t(coef(mlm_fit1_log_yolk)[4,]) %*% coef(mlm_fit1_log_yolk)[4,])

#code for coefficient of determination:
sum(diag(cov(Morph_log_yolk[,1:3])))
sum(diag(cov(mlm_fit1_log_yolk$fitted)))
sum(diag(cov(mlm_fit1_log_yolk$fitted)))/sum(diag(cov(Morph_log_yolk[,1:3])))
#model accounts for 57% of variance? seems high

#figure out if we need to do permutation test stuff to assess whether data conform to assumptions
#visualization:
dwplot(mlm_fit1_log_yolk) #this one doesn't work for some reason
plot(allEffects(mlm_fit1_log_yolk)) #this sort of works - maybe try to fix it up a bit
plot(emmeans(mlm_fit1_log_yolk, ~Treatment)) #this is useless
#is a ggplot possible?
 
#geomorph model:
mlm_fit2_log_yolk <- procD.lm(f1 = Morph_log_yolk[, 1:3] ~ Treatment*age, 
                         data = Morph_log_yolk, iter = 2000 )
summary(mlm_fit2_log_yolk)
coef(mlm_fit2_log_yolk)
#this basically gives same answer as first model

#create coefficient plots?