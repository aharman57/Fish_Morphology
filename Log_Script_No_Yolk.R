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
plot.(mlm_fit1_log)
summary(manova(mlm_fit1_log), test = "Wilks")
coef(mlm_fit1_log)
#would need to back-transform effect sizes to get to biologically relevant scale - see class lecture slides

#magnitude of treatment and age constrast vectors - but what does this really mean?
sqrt(t(coef(mlm_fit1_log)[2,]) %*% coef(mlm_fit1_log)[2,])
sqrt(t(coef(mlm_fit1_log)[3,]) %*% coef(mlm_fit1_log)[3,])
sqrt(t(coef(mlm_fit1_log)[4,]) %*% coef(mlm_fit1_log)[4,])

#code for coefficient of determination:
sum(diag(cov(Morph_log[,1:7])))
sum(diag(cov(mlm_fit1_log$fitted)))
sum(diag(cov(mlm_fit1_log$fitted)))/sum(diag(cov(Morph_log[,1:7])))
#model accounts for 66% of variance? seems high

#figure out if we need to do permutation test stuff to assess whether data conform to assumptions
#visualization:
dwplot(mlm_fit1_log) #this one doesn't work for some reason
plot(allEffects(mlm_fit1_log)) #this sort of works - maybe try to fix it up a bit - issue with ages
plot(emmeans(mlm_fit1_log,~Treatment)) #this is useless
#is a ggplot possible?

#geomorph model:
mlm_fit2_log <- procD.lm(f1 = Morph_log[, 1:7] ~ Treatment*age, 
                     data = Morph_log, iter = 2000 )
summary(mlm_fit2_log)
coef(mlm_fit2_log)
#this basically gives same answer as first model

#create coefficient plots?

#lmer model:
Morph_melt <- (Morph_clean_body
               %>% mutate(units=factor(1:n()))
               %>% gather(trait, value, -c(units, age, Treatment))
               %>% drop_na() #may not need this if we already omitted
               %>% arrange(units)
)

t1 <- system.time(
  lmer1 <- lmer(value ~ trait:(age*Treatment) - 1 +
                  (trait-1|units),
                data=Morph_melt,
                control=lmerControl(optCtrl=list(ftol_abs=1e-10),
                                    optimizer="bobyqa",
                                    check.nobs.vs.nlev="ignore",
                                    check.nobs.vs.nRE="ignore"))
)
summary(lmer1)
cc1 <- tidy(lmer1,effect="fixed") %>%
  tidyr::separate(term,into=c("trait","fixeff"),extra="merge",
                  remove=FALSE)
dwplot(cc1)+
  geom_vline(xintercept=0,lty=2) #this works but everything on different scales and fin values are super high
