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
          %>% select(-c(1:2, 16:17, 20:25)) #got rid of other variables we probably won't use
          %>% rename(Length = "Length (mm)",
                     Eye = "Measure Eye Diameter",
                     Fin_Anterior = "Measure Dorsal Fin Anterior Maximum",
                     Fin_Min = "Measure Dorsal Fin Minimum",
                     Fin_Posterior = "Measure Dorsal Fin Posterior Maximum",
                     Fin_Indent = "Fin indentation ratio",
                     Yolk_Width ="Measure Yolk Width",
                     Yolk_Height = "Measure Yolk Height",
                     Yolk_Vol = "Yolk volume (mm2)",
                     Jaw = "Jaw gape (um)",
                     Body_Weight = "body weight",
                     Yolk_Weight = "yolk weight"
          )
          %>% na.omit()
          %>% as.numeric(Fin_Anterior)
)

Morph_clean <- (Morph_clean
                %>% as.numeric("Fin_Anterior")
)
for (i in 1:13) {
  Morph_clean[,i] <- as.numeric(Morph_clean[,i])
}
as.numeric(Morph_clean$Body_Weight)
#covariance matrix of response traits:
cov(Morph[,3:11])
cor(Morph[,3:11])

pairs(Morph[, 3:11],
      pch = ".", gap = 0)

#scale response variables of interest:
Morphdata2 <- Morph %>%
  mutate(Length = scale(Length), Eye = scale(Eye), `Yolk volume` = scale(`Yolk volume`), Fin = scale(Fin), `Body weight` = scale(`Body weight`), Jaw = scale(Jaw), Yolk = scale(Yolk))

#drop variables we're not using (total weight and condition):
Morphdata2 <- select(Morphdata2, -c("Total weight", "Condition"))

Morph_melt <- (Morphdata2
             %>% mutate(units=factor(1:n()))
             %>% gather(trait,value, -c(units, Age, Treatment))
             %>% drop_na()
             %>% arrange(units)
)
#got this Warning message:
#attributes are not identical across measure variables;
#they will be dropped
#not sure if this is a problem?

#fit linear model:
t1 <- system.time(
  lmer1 <- lmer(value ~ trait:(Age*Treatment) - 1 +
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
