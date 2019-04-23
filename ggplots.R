library(tidyverse)
library(ggplot2)
load("Morph_clean.R")
load("Morph_log.R")
load("Morph_yolk.R")

#Looking at effect of age and treatment on body morphology and yolk variables:
length_boxplot <- (ggplot(Morph_clean_body, aes(x=age, y=Length, colour=Treatment))
                +geom_boxplot())
print(length_boxplot)

weight_boxplot <- (ggplot(Morph_clean_body, aes(x=age, y=Body_Weight, colour=Treatment))
                   +geom_boxplot())
print(weight_boxplot)

Fin_Ant_boxplot <- (ggplot(Morph_clean_body, aes(x=age, y=Fin_Anterior, colour=Treatment))
                   +geom_boxplot())
print(Fin_Ant_boxplot)

Fin_Min_boxplot <- (ggplot(Morph_clean_body, aes(x=age, y=Fin_Min, colour=Treatment))
                   +geom_boxplot())
print(Fin_Min_boxplot)

Fin_Post_boxplot <- (ggplot(Morph_clean_body, aes(x=age, y=Fin_Posterior, colour=Treatment))
                   +geom_boxplot())
print(Fin_post_boxplot)

Eye_boxplot <- (ggplot(Morph_clean_body, aes(x=age, y=Eye, colour=Treatment))
                   +geom_boxplot())
print(Eye_boxplot)

Jaw_boxplot <- (ggplot(Morph_clean_body, aes(x=age, y=Jaw, colour=Treatment))
                   +geom_boxplot())
print(Jaw_boxplot)
##change axis titles

yolkweight_boxplot <- (ggplot(Morph_clean_yolk, aes(x=age, y=Yolk_Weight, colour=Treatment))
                        +geom_boxplot())
print(yolkweight_boxplot)

yolkheight_boxplot <- (ggplot(Morph_clean_yolk, aes(x=age, y=Yolk_Height, colour=Treatment))
                       +geom_boxplot())
print(yolkheight_boxplot)

yolkwidth_boxplot <- (ggplot(Morph_clean_yolk, aes(x=age, y=Yolk_Width, colour=Treatment))
                       +geom_boxplot())
print(yolkwidth_boxplot)


#investigating relationships between length and other body morphology traits:

length_weight <- (ggplot(Morph_log, aes(x=Length, y=Body_Weight, colour=Treatment))
                     +geom_point()
                     +facet_grid(.~Treatment)
                     +geom_smooth(method="lm"))
print(length_weight)

length_FinAnt <- (ggplot(Morph_log, aes(x=Length, y=Fin_Anterior, colour=Treatment))
                     +geom_point()
                     +geom_smooth(method="lm")
                     +facet_grid(.~Treatment))
print(length_FinAnt)

length_FinMin <- (ggplot(Morph_log, aes(x=Length, y=Fin_Min, colour=Treatment))
                  +geom_point()
                  +geom_smooth(method="lm")
                  +facet_grid(.~Treatment))
print(length_FinMin)

length_FinPost <- (ggplot(Morph_log, aes(x=Length, y=Fin_Posterior, colour=Treatment))
                  +geom_point()
                  +geom_smooth(method="lm")
                  +facet_grid(.~Treatment))
print(length_FinPost)

length_eye <- (ggplot(Morph_log, aes(x=Length, y=Eye, colour=Treatment))
                  +geom_point()
                  +geom_smooth(method="lm")
                  +facet_grid(.~Treatment))
print(length_eye)

length_jaw <- (ggplot(Morph_log, aes(x=Length, y=Jaw, colour=Treatment))
                  +geom_point()
                  +geom_smooth(method="lm")
                  +facet_grid(.~Treatment))
print(length_jaw)

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
