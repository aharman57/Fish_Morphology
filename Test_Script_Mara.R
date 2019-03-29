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
          %>% select(-c(1:2, 4, 16:17, 20:25)) #got rid of other variables we probably won't use
          %>% rename(Length = "Length (mm)",
                     Eye = "Eye size (mm)",
                     Fin_Anterior = "Measure Dorsal Fin Anterior Maximum",
                     Fin_Min = "Measure Dorsal Fin Minimum",
                     Fin_Posterior = "Measure Dorsal Fin Posterior Maximum",
                     Fin_Indent = "Fin indentation ratio",
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

#next: need to convert treatment groups to factor or numeric after dealing with different names for same treatment


#covariance matrix of response traits:
cov(Morph_clean[,1:12])
cor(Morph_clean[,1:12])

pairs(Morph_clean[, 1:12],
      pch = ".", gap = 0)

#scale response variables of interest:
Morph_scaled <- (Morph_clean
                 %>% mutate(Length = scale(Length),
                            Eye = scale(Eye),
                            Fin_Anterior = scale(Fin_Anterior),
                            Fin_Min = scale(Fin_Min),
                            Fin_Posterior = scale(Fin_Posterior),
                            Fin_Indent = scale(Fin_Indent),
                            Yolk_Width = scale(Yolk_Width),
                            Yolk_Height = scale(Yolk_Height),
                            Yolk_Vol = scale(Yolk_Vol),
                            Jaw = scale(Jaw),
                            Body_Weight = scale(Body_Weight),
                            Yolk_Weight = scale(Yolk_Weight)
                 )
)

mlm_fit1 <- lm(as.matrix(Morph_scaled[,1:12]) ~ Treatment*age, data = Morph_scaled)
class(mlm_fit1) #not sure if we need this
summary(manova(mlm_fit1), test = "Wilks")



#for lmer (but we may just be able to use simple linear model from first lecture slide?)
Morph_melt <- (Morph_scaled
             %>% mutate(units=factor(1:n()))
             %>% gather(trait,value, -c(units, age, Treatment))
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

### is this relevant/important??? ###
covariance1 = cov(Morph)
print(covariance1)

### correlation matrix for all predictor variables ###
correlation1 = cor(Morph)
print(correlation1)
### shows that body weight and eye diameter are most correlated variables ###