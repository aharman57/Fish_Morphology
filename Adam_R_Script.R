library(tidyverse)
library(MCMCglmm)
library(lme4)
library(corrplot)
library(broom.mixed)
library(dotwhisker)
library(ggplot2); theme_set(theme_bw())

Morph <- (read_csv(file="Morph_Data_2016-2017.csv"))
Morph <- select(Morph,-c(20:25))

Morph_clean <- na.omit(Morph)

Morph_clean <- as.numeric(Morph_clean$`Measure Dorsal Fin Anterior Maximum`)
Morph_clean <- as.numeric(Morph_clean$`Measure Dorsal Fin Minimum`)
Morph_clean <- as.numeric(Morph_clean$`Measure Dorsal Fin Posterior Maximum`)
Morph_clean <- as.numeric(Morph_clean$`Measure Yolk Width`)
Morph_clean <- as.numeric(Morph_clean$`Measure Yolk Height`)
       
### is this relevant/important??? ###
covariance1 = cov(Morph[,3:17])
print(covariance1)

### correlation matrix for all predictor variables ###
correlation1 = cor(Morph)
print(correlation1)
### shows that body weight and eye diameter are most correlated variables ###

pairs(Morph,
      pch = ".", gap = 0)

