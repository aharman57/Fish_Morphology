#combining yolk and body size into one model

library(tidyverse)
library(MCMCglmm)
library(lme4)
library(corrplot)
library(broom.mixed)
library(dotwhisker)
library(ggplot2); theme_set(theme_bw())
library(car)
library(geomorph)
library(emmeans)
library(dotwhisker)
library(effects)
library(MASS)

Morph <- read_csv("Morph_Data_2016-2017.csv")
Morph_clean_body_yolk <- (Morph
                     %>% select(-c(1:13, 16:17, 20:25))
                     %>% rename(Body_Weight = "body weight",
                                Yolk_Weight = "yolk weight",
                                Treatment = "treatment group"
                     )
                     %>% na.omit()
)

Morph_clean_body_yolk$Treatment <- as.factor(Morph_clean_body$Treatment)
Morph_clean_body_yolk$age <- as.factor(Morph_clean_body$age)
Morph_log_body_yolk <- (Morph_clean_body_yolk
                        %>% mutate(Body_Weight = log(Body_Weight),
                                       Yolk_Weight = log(Yolk_Weight)
                                   )
)