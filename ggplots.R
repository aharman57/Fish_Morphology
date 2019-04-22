library(tidyverse)
library(ggplot2)
load("Morph_clean.R")
load("Morph_log.R")
load("Morph_yolk.R")

length_boxplot <- (ggplot(Morph_log, aes(x=age, y=Length, colour=Treatment))
                +geom_boxplot())
print(length_boxplot)

yolkweight_boxplot <- (ggplot(Morph_clean_yolk, aes(x=age, y=Yolk_Weight, colour=Treatment))
                        +geom_boxplot())
print(yolkweight_boxplot)

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

#need to change dataset
#body_yolk <- (ggplot(Morph, aes(x="body weight", y="yolk weight", colour=age))
#              +geom_point()
#              +facet_grid(.~"treatment group"))
#print(body_yolk)
