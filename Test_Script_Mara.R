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
Morph_clean <- (Morph
          %>% select(-c(1:2, 4, 9, 16:17, 20:25)) #got rid of other variables we probably won't use
          %>% rename(Length = "Length (mm)",
                     Eye = "Eye size (mm)",
                     Fin_Anterior = "Measure Dorsal Fin Anterior Maximum",
                     Fin_Min = "Measure Dorsal Fin Minimum",
                     Fin_Posterior = "Measure Dorsal Fin Posterior Maximum",
                     Yolk_Width ="Measure Yolk Width",
                     Yolk_Height = "Measure Yolk Height",
                     Yolk_Vol = "Yolk volume (mm2)",
                     Jaw = "Jaw gape (um)",
                     Body_Weight = "body weight",
                     Yolk_Weight = "yolk weight",
                     Treatment = "treatment group"
          )
          %>% na.omit()
)

Morph_clean$Fin_Anterior <- as.numeric(Morph_clean$Fin_Anterior)
Morph_clean$Fin_Min <- as.numeric(Morph_clean$Fin_Min)
Morph_clean$Fin_Posterior <- as.numeric(Morph_clean$Fin_Posterior)
Morph_clean$Yolk_Width <- as.numeric(Morph_clean$Yolk_Width)
Morph_clean$Yolk_Height <- as.numeric(Morph_clean$Yolk_Height)


#covariance matrix of response traits:
cov(Morph_clean[,1:11])
cor(Morph_clean[,1:11])

pairs(Morph_clean[, 1:11],
      pch = ".", gap = 0)

#this plot is weird and not working - maybe remove
scatterplotMatrix( ~ Length + Eye + Fin_Anterior + Fin_Min | interaction(Treatment, age), 
                   ellipse = TRUE, data = Morph_clean, gap = 0,
                   plot.points = T, pch = 20, cex  = 0.5)

#should we have scaled first?
eig_vals <- svd(cov(Morph_clean[, 1:11]))$d
prod(eig_vals)
det(cov(Morph_clean[, 1:11]))
sum(eig_vals)
sum(diag(cov(Morph_clean[, 1:11])))

#scale response variables of interest:
Morph_scaled <- (Morph_clean
                 %>% mutate(Length = scale(Length),
                            Eye = scale(Eye),
                            Fin_Anterior = scale(Fin_Anterior),
                            Fin_Min = scale(Fin_Min),
                            Fin_Posterior = scale(Fin_Posterior),
                            Yolk_Width = scale(Yolk_Width),
                            Yolk_Height = scale(Yolk_Height),
                            Yolk_Vol = scale(Yolk_Vol),
                            Jaw = scale(Jaw),
                            Body_Weight = scale(Body_Weight),
                            Yolk_Weight = scale(Yolk_Weight)
                 )
)

#eigenvalues after scaling:
eig_vals_scaled <- svd(cov(Morph_scaled[, 1:11]))$d
prod(eig_vals_scaled)
det(cov(Morph_scaled[, 1:11]))
sum(eig_vals_scaled)
sum(diag(cov(Morph_scaled[, 1:11])))

mlm_fit1 <- lm(as.matrix(Morph_scaled[,1:11]) ~ Treatment*age, data = Morph_scaled)
summary(manova(mlm_fit1), test = "Wilks")
coef(mlm_fit1)

#magnitude of treatment and age constrast vectors - but what does this really mean?
sqrt(t(coef(mlm_fit1)[2,]) %*% coef(mlm_fit1)[2,])
sqrt(t(coef(mlm_fit1)[3,]) %*% coef(mlm_fit1)[3,])
sqrt(t(coef(mlm_fit1)[4,]) %*% coef(mlm_fit1)[4,])

#code for coefficient of determination:
sum(diag(cov(Morph_scaled[,1:11])))
sum(diag(cov(mlm_fit1$fitted)))
sum(diag(cov(mlm_fit1$fitted)))/sum(diag(cov(Morph_scaled[,1:11])))
#seems like a very high number

#geomorph model:
mlm_fit2 <- procD.lm(f1 = Morph_scaled[, 1:11] ~ Treatment*age, 
                     data = Morph_scaled, iter = 2000 )
summary(mlm_fit2)
coef(mlm_fit2)

#also need to check model assumptions somehow?

#for lmer (but we may just be able to use simple linear model from first lecture slide?)
Morph_melt <- (Morph_scaled
             %>% mutate(units=factor(1:n()))
             %>% gather(trait, value, -c(units, age, Treatment))
             %>% drop_na() #may not need this if we already omitted
             %>% arrange(units)
)
#got this Warning message:
#attributes are not identical across measure variables;
#they will be dropped
#not sure if this is a problem?

#fit linear model (this didn't work):
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
