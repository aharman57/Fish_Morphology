## Fish_Morphology: Final Project for Biology 708

### Data

Our dataset consists of information on several morphological measurements of Lake Whitefish hatchlings that were reared in different water temperatures. The Age represents the number of days after hatching at which morphological measurements were taken. Treatment values represent the temperature (in degrees Celsius) of the water in which eggs were incubated. After the incubation treatment (post-hatch) larvae were all kept at the same temperature (8 degC). Spawning adults were caught in Lake Huron and stripped of eggs/milt, and fertilized embryos brought back to the lab. Sampling was lethal, so each row represents a different fish (different fish were sampled at different ages). This experiment was repeated several years in a row. This is only 1 year of the dataset, so we are working on acquiring several more years of data. All the response variables are correlated, as they all come from the same individual. As the fish grows, so do all the other body features that we measure (except the yolk, which diminishes). A simple linear model won't account for the fact that these variables are correlated, so I want to be able to fit a model that accounts for the correlations.

### Scientific Questions

1. What is the influence of temperature on length and biomass? How does temperature influence the relationship between length and biomass?

2. What is the relationship between length and other morphological characteristics, and how does temperature affect these relationships?

3. What is the effect of temperature on 'yolk efficiency' (i.e. relationship between growth rate (increase in length or biomass) and reduction in yolk volume/mass over time)?

*BMB: good for you for not asking "is there a difference due to ... " ?  Could you make more precise predictions (e.g. do you expect yolk efficiency to decrease with temperature?)*

### Analysis Plans

- Calculations of growth rate and yolk depletion rate over different time periods for each treatment and year to be used in analysis of yolk efficiency
- Do correlation tests to see how characters are correlated *BMB: this sounds exploratory (although interesting): you don't need to know in advance how they're correlated to do the next step.  See the `corrplot` package for pretty correlation pictures (also `pairs(X, pch=".", gap=0)` or use the scatterplot/pairs function that Ian Dworkin showed in his lecture)*
- Come up with a model (potentially a mixed linear model) to test effects of temperature on morphological characteristics, accounting for correlation between variables

Questions/Issues:

1. How do we deal with correlation between different morphological characteristics? Correlation makes it difficult to tell whether temperature is directly affecting a certain morphological feature, or if this is just due to correlation with another morphological feature.

*BMB: there's basically no way to pull out causality from a statistical analysis (you can use structural equation models, but this basically means setting some alternative hypotheses about how variables are correlated - you're not actually getting causal inference for nothing)*

2. How do we test interactions between multiple morphological characteristics so that we can understand if temperature influences the relationships between these characters?

*BMB: interesting question. Most multivariate models assume a **constant** correlation structure. The easiest thing is to see how the temperature coefficients for different traits vary; if trait parameters differ, then the allometry changes with temperature [in principle you can even ask whether they're sig diff from each other ...]*

3. Can we create a model that describes the relationships between different morphological characteristics so that we can infer the response of one character based on the data about other characters? We want to be able to accurately predict total length based on a few morphological measurements, to make data collection easier.

*BMB: also interesting. Given a fitted multivariate model you should be able to get an estimate of one trait from a subset of the other traits, in principle (I'll have to think a little bit more about how to actually do it: in spatial modeling, this is called **kriging**)*
