### work in Adam branch, type in terminal - git merge Adam ###


library(tidyverse)
library(MCMCglmm)
library(lme4)
library(corrplot)
library(broom.mixed)
library(dotwhisker)
library(ggplot2); theme_set(theme_bw())

Morph <- (read_csv(file="Morph_Data.csv")
          %>% rename(Age="Age (dph)",
                     Eye="Eye diameter",
                     Fin="Fin indent",
                     Yolk="Yolk weight",
                     Jaw="Jaw Gape")
)

Treatment <- Morph$Treatment
Age <- Morph$Age
Eye <- Morph$Eye
Fin <- Morph$Fin
Yolk <- Morph$Yolk
Jaw <- Morph$Jaw
Length <- Morph$Length


### is this relevant/important??? ###
covariance1 = cov(Morph)
print(covariance1)

### correlation matrix for all predictor variables ###
correlation1 = cor(Morph)
print(correlation1)
### shows that body weight and eye diameter are most correlated variables ###
