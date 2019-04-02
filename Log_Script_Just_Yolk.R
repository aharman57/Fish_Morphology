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
                     %>% select(-c(1:9, 13:14, 16:17, 20:25)) #got rid of other variables we probably won't use
                     %>% rename(Yolk_Width ="Measure Yolk Width",
                                Yolk_Height = "Measure Yolk Height",
                                Yolk_Vol = "Yolk volume (mm2)",
                                Yolk_Weight = "yolk weight",
                                Treatment = "treatment group"
                     )
                     %>% na.omit()
                     %>% filter(age != 28) #may need to filter out 14 days too or else remove yolk weight
)

#covariance matrix of response traits:
cov(Morph_clean_yolk[,1:4])
cor(Morph_clean_yolk[,1:4])

pairs(Morph_clean_yolk[, 1:4],
      pch = ".", gap = 0)

#logging variables instead of scaling:
Morph_log_yolk <- (Morph_clean_yolk
              %>% mutate(Yolk_Width = log(Yolk_Width),
                         Yolk_Height = log(Yolk_Height),
                         Yolk_Vol = log(Yolk_Vol),
                         Yolk_Weight = log(Yolk_Weight)
              )
)

pairs(Morph_log_yolk[, 1:4],
      pch = ".", gap = 0)

#eigenvalues after logging - need to deal with zeros in weight first:
eig_vals_log_yolk <- svd(cov(Morph_log_yolk[, 1:4]))$d
prod(eig_vals_log_yolk)
sum(eig_vals_log_yolk)

mlm_fit1_log_yolk <- lm(as.matrix(Morph_log_yolk[,1:4]) ~ Treatment*age, data = Morph_log_yolk)
summary(manova(mlm_fit1_log_yolk), test = "Wilks")
coef(mlm_fit1_log_yolk)
#would need to back-transform effect sizes to get to biologically relevant scale - see class lecture slides

#magnitude of treatment and age constrast vectors - but what does this really mean?
sqrt(t(coef(mlm_fit1_log_yolk)[2,]) %*% coef(mlm_fit1_log_yolk)[2,])
sqrt(t(coef(mlm_fit1_log_yolk)[3,]) %*% coef(mlm_fit1_log_yolk)[3,])
sqrt(t(coef(mlm_fit1_log_yolk)[4,]) %*% coef(mlm_fit1_log_yolk)[4,])

#code for coefficient of determination:
sum(diag(cov(Morph_log_yolk[,1:4])))
sum(diag(cov(mlm_fit1_log_yolk$fitted)))
sum(diag(cov(mlm_fit1_log_yolk$fitted)))/sum(diag(cov(Morph_log_yolk[,1:4])))
#model accounts for 66% of variance? seems high

#figure out if we need to do permutation test stuff to assess whether data conform to assumptions

#geomorph model:
mlm_fit2_log_yolk <- procD.lm(f1 = Morph_log_yolk[, 1:4] ~ Treatment*age, 
                         data = Morph_log_yolk, iter = 2000 )
summary(mlm_fit2_log_yolk)
coef(mlm_fit2_log_yolk)
#this basically gives same answer as first model

#create coefficient plots?