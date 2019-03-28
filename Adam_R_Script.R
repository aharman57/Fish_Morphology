library(tidyverse)
library(MCMCglmm)
library(lme4)
library(corrplot)
library(broom.mixed)
library(dotwhisker)
library(ggplot2); theme_set(theme_bw())

Morph <- (read_csv(file="Morph_Data_2016-2017.csv"))
         
### is this relevant/important??? ###
covariance1 = cov(Morph)
print(covariance1)

### correlation matrix for all predictor variables ###
correlation1 = cor(Morph)
print(correlation1)
### shows that body weight and eye diameter are most correlated variables ###

pairs(Morph,
      pch = ".", gap = 0)

