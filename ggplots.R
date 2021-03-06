library(tidyverse)
library(ggplot2)
load("Morph_clean.R")
load("Morph_log.R")
load("Morph_yolk.R")

#Looking at effect of age and treatment on body morphology and yolk variables:
length_boxplot <- (ggplot(Morph_clean_body, aes(x=age, y=Length, colour=Treatment))
                +geom_boxplot()
                +labs(x="Age (days post-hatch)", y="Body Length (mm)", colour= "Incubation Temperature \n(degrees C)"))
print(length_boxplot)

weight_boxplot <- (ggplot(Morph_clean_body, aes(x=age, y=Body_Weight, colour=Treatment))
                   +geom_boxplot()
                   +labs(x="Age (days post-hatch)", y="Body Weight (g)", colour="Incubation Temperature \n(degrees C)"))
print(weight_boxplot)

Fin_Ant_boxplot <- (ggplot(Morph_clean_body, aes(x=age, y=Fin_Anterior, colour=Treatment))
                   +geom_boxplot()
                   +labs(x="Age (days post-hatch)", y="Anterior Fin Width (mm)", colour= "Incubation Temperature \n(degrees C)"))
print(Fin_Ant_boxplot)

Fin_Min_boxplot <- (ggplot(Morph_clean_body, aes(x=age, y=Fin_Min, colour=Treatment))
                   +geom_boxplot()
                   +labs(x="Age (days post-hatch)", y="Minimum Fin Width (mm)", colour= "Incubation Temperature \n(degrees C)"))
print(Fin_Min_boxplot)

Fin_Post_boxplot <- (ggplot(Morph_clean_body, aes(x=age, y=Fin_Posterior, colour=Treatment))
                   +geom_boxplot()
                   +labs(x="Age (days post-hatch)", y="Posterior Fin Width (mm)", colour= "Incubation Temperature \n(degrees C)"))
print(Fin_Post_boxplot)

Eye_boxplot <- (ggplot(Morph_clean_body, aes(x=age, y=Eye, colour=Treatment))
                   +geom_boxplot()
                +labs(x="Age (days post-hatch)", y="Eye Diameter (mm)", colour= "Incubation Temperature \n(degrees C)"))
print(Eye_boxplot)

Jaw_boxplot <- (ggplot(Morph_clean_body, aes(x=age, y=Jaw, colour=Treatment))
                   +geom_boxplot()
                +labs(x="Age (days post-hatch)", y="Jaw Gape (um)", colour= "Incubation Temperature \n(degrees C)"))
print(Jaw_boxplot)

yolkweight_boxplot <- (ggplot(Morph_clean_yolk, aes(x=age, y=Yolk_Weight, colour=Treatment))
                        +geom_boxplot()
                       +labs(x="Age (days post-hatch)", y="Yolk Weight (g)", colour= "Incubation Temperature \n(degrees C)"))
print(yolkweight_boxplot)

yolkheight_boxplot <- (ggplot(Morph_clean_yolk, aes(x=age, y=Yolk_Height, colour=Treatment))
                       +geom_boxplot()
                       +labs(x="Age (days post-hatch)", y="Yolk Height (mm)", colour= "Incubation Temperature \n(degrees C)"))
print(yolkheight_boxplot)

yolkwidth_boxplot <- (ggplot(Morph_clean_yolk, aes(x=age, y=Yolk_Width, colour=Treatment))
                       +geom_boxplot()
                      +labs(x="Age (days post-hatch)", y="Yolk Width (mm)", colour= "Incubation Temperature \n(degrees C)"))
print(yolkwidth_boxplot)


#investigating relationships between length and other body morphology traits:

length_weight <- (ggplot(Morph_log, aes(x=Length, y=Body_Weight, colour=Treatment))
                     +geom_point()
                     +geom_smooth(method="lm")
                    +facet_grid(.~Treatment)
                    +labs(x="Body Length (mm)", y="Body Weight (g)", colour= "Incubation Temperature \n(degrees C)"))
print(length_weight)

length_FinAnt <- (ggplot(Morph_log, aes(x=Length, y=Fin_Anterior, colour=Treatment))
                     +geom_point()
                     +geom_smooth(method="lm")
                     +facet_grid(.~Treatment)
                  +labs(x="Body Length (mm)", y="Anterior Fin Width (mm)", colour= "Incubation Temperature \n(degrees C)"))
print(length_FinAnt)

length_FinMin <- (ggplot(Morph_log, aes(x=Length, y=Fin_Min, colour=Treatment))
                  +geom_point()
                  +geom_smooth(method="lm")
                  +facet_grid(.~Treatment)
                  +labs(x="Body Length", y="Minimum Fin Width (mm)", colour= "Incubation Temperature \n(degrees C)"))
print(length_FinMin)

length_FinPost <- (ggplot(Morph_log, aes(x=Length, y=Fin_Posterior, colour=Treatment))
                  +geom_point()
                  +geom_smooth(method="lm")
                  +facet_grid(.~Treatment)
                  +labs(x="Body Length (mm)", y="Posterior Fin Width (mm)", colour= "Incubation Temperature \n(degrees C)"))
print(length_FinPost)

length_eye <- (ggplot(Morph_log, aes(x=Length, y=Eye, colour=Treatment))
                  +geom_point()
                  +geom_smooth(method="lm")
                  +facet_grid(.~Treatment)
               +labs(x="Body Length (mm)", y="Eye Diameter (mm)", colour= "Incubation Temperature \n(degrees C)"))
print(length_eye)

length_jaw <- (ggplot(Morph_log, aes(x=Length, y=Jaw, colour=Treatment))
                  +geom_point()
                  +geom_smooth(method="lm")
                  +facet_grid(.~Treatment)
                  +labs(x="Body Length (mm)", y="Jaw Gape (um)", colour= "Incubation Temperature \n(degrees C)"))
print(length_jaw)

#Investigating relationship between body weight and yolk across treatments and ages
load("Morph_bodyyolk")
body_yolk <- (ggplot(Morph_clean_body_yolk, aes(x=Body_Weight, y=Yolk_Weight, colour=age))
             +geom_point()
             +facet_grid(.~Treatment)
             +geom_smooth(method="lm")
             +labs(x="Body Weight (g)", y="Yolk Weight (g)", colour= "Age \n(days post-hatch)"))
print(body_yolk)
body_yolk2 <- (ggplot(Morph_clean_body_yolk, aes(x=Body_Weight, y=Yolk_Weight, colour=Treatment))
               +geom_point(size = 0.5)
               +facet_grid(.~age, labeller = label_both)
               +geom_smooth(method="lm")
               +labs(x="Body Weight (g)", y="Yolk Weight (g)", colour= "Incubation Temperature \n(degrees C)"))
print(body_yolk2)
