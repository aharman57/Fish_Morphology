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

morph_vars <- c("Length","tibia","tarsus","SCT")
morph_vars_sc <- paste(morph_vars,"s",sep="_")
dll_data2 <- dll_data
### is this relevant/important??? ###
covariance1 = cov(Morph)
print(covariance1)

### correlation matrix for all predictor variables ###
correlation1 = cor(Morph)
print(correlation1)
### shows that body weight and eye diameter are most correlated variables ###
