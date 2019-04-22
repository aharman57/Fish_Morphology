library(tidyverse)
library(ggplot2)
load("Morph_clean.R")
load("Morph_log.R")
load("Morph_yolk.R")

#Looking at effect of age and treatment on body morphology and yolk variables:
length_boxplot <- (ggplot(Morph_log, aes(x=age, y=Length, colour=Treatment))
                +geom_boxplot())
print(length_boxplot)
##do for each response, change axis titles, use unlogged data

yolkweight_boxplot <- (ggplot(Morph_clean_yolk, aes(x=age, y=Yolk_Weight, colour=Treatment))
                        +geom_boxplot())
print(yolkweight_boxplot)


#investigating relationships between length and other body morphology traits:
length_vs_FinAnt <- (ggplot(Morph_log, aes(x=Length, y=Fin_Anterior, colour=Treatment))
                     +geom_point()
                     +geom_smooth(method="lm")
                     +facet_grid(.~Treatment))
print(length_vs_FinAnt)
length_vs_weight <- (ggplot(Morph_log, aes(x=Length, y=Body_Weight, colour=Treatment))
                     +geom_point()
                     +facet_grid(.~Treatment)
                     +geom_smooth(method="lm"))
print(length_vs_weight)
#do all and see which ones are most different

#Investigating relationship between body weight and yolk across treatments and ages
load("Morph_bodyyolk")
body_yolk <- (ggplot(Morph_clean_body_yolk, aes(x=Body_Weight, y=Yolk_Weight, colour=age))
             +geom_point()
             +facet_grid(.~Treatment)
             +geom_smooth(method="lm", se=FALSE))
print(body_yolk)
body_yolk2 <- (ggplot(Morph_clean_body_yolk, aes(x=Body_Weight, y=Yolk_Weight, colour=Treatment))
               +geom_point()
               +facet_grid(.~age)
               +geom_smooth(method="lm", se=FALSE))
print(body_yolk2)