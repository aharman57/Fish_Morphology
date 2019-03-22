library(tidyverse)

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
covariance1 = cov(Length,Age)
print(covariance1)
covariance2 = cov(Length,Treatment)
print(covariance2)
covariance3 = cov(Length,Yolk)
print(covariance3)

### correlation matrix for all predictor variables ###
correlation1 = cor(Morph)
print(correlation1)
### shows that body weight and eye diameter are most correlated variables ###
