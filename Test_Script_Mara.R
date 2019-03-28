library(tidyverse)
library(MCMCglmm)
library(lme4)
library(corrplot)
library(broom.mixed)
library(dotwhisker)
library(ggplot2); theme_set(theme_bw())
library(car)
library(geomorph)

Morph <- (read_csv(file="Morph_Data.csv")
          %>% rename(Age="Age (dph)",
                     Eye="Eye diameter",
                     Fin="Fin indent",
                     Yolk="Yolk weight",
                     Jaw="Jaw Gape")
)

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
